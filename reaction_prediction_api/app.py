from flask import Flask, render_template, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AllChem import ReactionFromSmarts
import json

app = Flask(__name__)

# Your reaction library
reaction_library = {
    'amide_formation': '[C:1](=[O:2])-[OD1:3].[N!H0:4]>>[C:1](=[O:2])[N:4]',
    'esterification': '[C:1](=[O:2])-[OD1:3].[O!H0:4]-[C:5]>>[C:1](=[O:2])[O:4][C:5]',
    'nitration': '[cH:1]>>[c:1][N+](=O)[O-]',
    'hydrogenation': '[C:1]=[C:2]>>[C:1][C:2]',
    'ester_hydrolysis': '[O:1]=[C:2][O:3][C:4]>>[O:1]=[C:2][OH].[OH:3][C:4]',
    'reductive_amination': '[C:1]=[O:2].[N!H0:3]>>[C:1][N:3]',
}

def run_reaction(reaction_smarts, reactant_list):
    """Applies a reaction and returns unique products."""
    rxn = ReactionFromSmarts(reaction_smarts)
    products_list = rxn.RunReactants(reactant_list)
    unique_products = []
    seen_smiles = set()
    
    for products in products_list:
        for prod in products:
            try:
                Chem.SanitizeMol(prod)
                prod_smiles = Chem.MolToSmiles(prod, canonical=True)
                if prod_smiles not in seen_smiles:
                    seen_smiles.add(prod_smiles)
                    unique_products.append(prod)
            except:
                continue
    return unique_products

def advanced_evaluation(reactants, products, reaction_name):
    """A more sophisticated scoring function."""
    if not products:
        return 0.0, ["No reaction"]
    
    score = 0.0
    reasons = []
    
    # Reaction-specific scoring
    if reaction_name == 'amide_formation':
        has_acid = any('[OD1]' in Chem.MolToSmarts(mol) for mol in reactants)
        has_amine = any(mol.HasSubstructMatch(Chem.MolFromSmarts('[N!H0]')) for mol in reactants)
        if has_acid and has_amine:
            score += 2.0
            reasons.append("Appropriate reactants for amide formation")
    
    elif reaction_name == 'esterification':
        has_acid = any('[OD1]' in Chem.MolToSmarts(mol) for mol in reactants)
        has_alcohol = any(mol.HasSubstructMatch(Chem.MolFromSmarts('[O!H0]')) for mol in reactants)
        if has_acid and has_alcohol:
            score += 2.0
            reasons.append("Appropriate reactants for esterification")
    
    # General scoring criteria
    product_mw = sum(Descriptors.MolWt(p) for p in products)
    reactant_mw = sum(Descriptors.MolWt(r) for r in reactants)
    
    # Mass balance check
    if reactant_mw > 0:  # Avoid division by zero
        mass_balance_ratio = product_mw / reactant_mw
        if 0.8 < mass_balance_ratio < 1.2:
            score += 1.0
            reasons.append("Good mass balance")
        else:
            score -= 1.0
            reasons.append(f"Poor mass balance: ratio {mass_balance_ratio:.2f}")
    
    return max(0.0, score), reasons

def find_best_reactions(reactant_smiles_list):
    """Try all reactions in the library and return scored results."""
    reactants = [Chem.MolFromSmiles(smi) for smi in reactant_smiles_list if smi.strip()]
    
    # Filter out None values (invalid SMILES)
    reactants = [mol for mol in reactants if mol is not None]
    
    if not reactants:
        return {"error": "No valid reactants provided"}
    
    results = []
    
    for rxn_name, rxn_smarts in reaction_library.items():
        try:
            products = run_reaction(rxn_smarts, reactants)
            score, reasons = advanced_evaluation(reactants, products, rxn_name)
            
            # Convert products to serializable format
            product_smiles = [Chem.MolToSmiles(prod) for prod in products]
            
            results.append({
                'reaction_name': rxn_name,
                'score': score,
                'products': product_smiles,
                'reasons': reasons,
                'smarts': rxn_smarts
            })
            
        except Exception as e:
            results.append({
                'reaction_name': rxn_name,
                'score': 0.0,
                'products': [],
                'reasons': [f"Error: {str(e)}"],
                'smarts': rxn_smarts
            })
    
    # Sort by score descending
    results.sort(key=lambda x: x['score'], reverse=True)
    return results

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict_reactions():
    try:
        data = request.get_json()
        reactant1 = data.get('reactant1', '').strip()
        reactant2 = data.get('reactant2', '').strip()
        
        # Create list of non-empty reactants
        reactants = [r for r in [reactant1, reactant2] if r]
        
        if not reactants:
            return jsonify({"error": "Please provide at least one reactant"})
        
        results = find_best_reactions(reactants)
        
        # Filter out reactions with score <= 0
        valid_results = [r for r in results if r['score'] > 0]
        
        return jsonify({
            "success": True,
            "reactants": reactants,
            "reactions": valid_results,
            "total_reactions_found": len(valid_results)
        })
        
    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"})

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
