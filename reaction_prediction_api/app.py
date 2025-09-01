from flask import Flask, render_template, request, jsonify, send_from_directory
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolToSmiles
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
import json

app = Flask(__name__)


# Serve static files from the static directory
@app.route('/static/<path:filename>')
def serve_static(filename):
    return send_from_directory('static', filename)

# Serve Kekule.js files specifically
@app.route('/kekule_libs/<path:filename>')
def serve_kekule(filename):
    return send_from_directory('static/kekule_libs', filename)



# Your reaction library
reaction_library = {
    'amide_formation': '[C:1](=[O:2])-[OD1:3].[N!H0:4]>>[C:1](=[O:2])[N:4]',
    'esterification': '[C:1](=[O:2])-[OD1:3].[O!H0:4]-[C:5]>>[C:1](=[O:2])[O:4][C:5]',
    'nitration': '[cH:1]>>[c:1][N+](=O)[O-]',
    'hydrogenation': '[C:1]=[C:2]>>[C:1][C:2]',
    'ester_hydrolysis': '[O:1]=[C:2][O:3][C:4]>>[O:1]=[C:2][OH].[OH:3][C:4]',
    'reductive_amination': '[C:1]=[O:2].[N!H0:3]>>[C:1][N:3]',
}



def smiles_to_cml_with_3d(smiles_string):
    """
    Generates a 3D-embedded molecule and returns its structure in CML format.
    Returns None if conversion fails.
    """
    try:
        if not smiles_string or not isinstance(smiles_string, str):
            print(f"Invalid SMILES input: {smiles_string}")
            return None
            
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print(f"RDKit failed to parse SMILES: {smiles_string}")
            return None

        # Add hydrogens for better 3D structure
        #mol_with_hs = Chem.AddHs(mol)
        
        AllChem.Compute2DCoords(mol)
        

        
        
        # Convert to CML format using RDKit's CML writer
        cml_block = Chem.MolToCMLBlock(mol)
        
        if cml_block.startswith('<?xml'):
            # Find the first line break after the XML declaration
            lines = cml_block.split('\n', 1)
            if len(lines) > 1:
                cml_block = lines[1]  # Keep everything after the first line
        
        print(cml_block)
        
        return cml_block
        
    except Exception as e:
        print(f"Error processing SMILES '{smiles_string}' with RDKit: {e}")
        return None

def smiles_to_cml_with_2d(smiles_string):
    """
    Generates a molecule with 2D coordinates and returns its structure in CML format.
    Returns None if conversion fails.
    """
    try:
        if not smiles_string or not isinstance(smiles_string, str):
            print(f"Invalid SMILES input: {smiles_string}")
            return None
            
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print(f"RDKit failed to parse SMILES: {smiles_string}")
            return None

        # Generate 2D coordinates (this is key for 2D display)
        AllChem.Compute2DCoords(mol)
        
        # Generate CML manually with 2D coordinates
        conf = mol.GetConformer()
        cml_lines = []
        cml_lines.append('<cml xmlns="http://www.xml-cml.org/schema">')
        cml_lines.append('  <molecule>')
        
        # Add atom array with 2D coordinates
        cml_lines.append('    <atomArray>')
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            element = atom.GetSymbol()
            # Use x2 and y2 for 2D coordinates instead of x3, y3, z3
            cml_lines.append(f'      <atom id="a{i}" elementType="{element}" x2="{pos.x:.6f}" y2="{pos.y:.6f}"/>')
        cml_lines.append('    </atomArray>')
        
        # Add bond array
        cml_lines.append('    <bondArray>')
        for i, bond in enumerate(mol.GetBonds()):
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            order = bond.GetBondTypeAsDouble()
            # Convert bond order to CML format
            if order == 1.0:
                bond_order = "S"
            elif order == 2.0:
                bond_order = "D"
            elif order == 3.0:
                bond_order = "T"
            elif order == 1.5:  # Aromatic
                bond_order = "A"
            else:
                bond_order = "S"  # Default to single
                
            cml_lines.append(f'      <bond id="b{i}" atomRefs2="a{begin_idx} a{end_idx}" order="{bond_order}"/>')
        cml_lines.append('    </bondArray>')
        
        cml_lines.append('  </molecule>')
        cml_lines.append('</cml>')
        
        return '\n'.join(cml_lines)
        
    except Exception as e:
        print(f"Error processing SMILES '{smiles_string}' with RDKit: {e}")
        return None


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
            

            
            product_cmls = []
            for smiles in product_smiles:
                cml = smiles_to_cml_with_2d(smiles)  # Passing individual SMILES string
                product_cmls.append(cml)
            
            
            
            results.append({
                'reaction_name': rxn_name,
                'score': score,
                'products': product_smiles,
                'products_cmls': product_cmls,
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

@app.route('/reactions', methods=['GET'])
def get_reactions():
    # Convert reaction_library to a list of dictionaries
    reactions = [{"id": key, "name": key.replace('_', ' ').capitalize()} for key in reaction_library.keys()]
    return jsonify(reactions)


@app.route('/predict', methods=['POST'])
def predict_reactions():
    try:
        data = request.get_json()
        print(data)
        reactant1 = data.get('reactant1', '').strip()
        reactant2 = data.get('reactant2', '').strip()
        selected_reaction = data.get('reaction', '').strip()  # Get selected reaction
        
        # Create list of non-empty reactants
        reactants = [r for r in [reactant1, reactant2] if r]
        
        if not reactants:
            return jsonify({"error": "Please provide at least one reactant"})
        
        # If a specific reaction is selected, run only that reaction
        print(f"Selected reaction: {selected_reaction}");
        if selected_reaction and selected_reaction in reaction_library:
            reaction_smarts = reaction_library[selected_reaction]
            print(f"Reaction Smarts: {reaction_smarts}")
            #products = run_reaction(rxn_smarts, reactants)
            #score, reasons = advanced_evaluation(reactants, products, rxn_name)
            print(f"Reactants: {reactants}")
            products = run_reaction(reaction_smarts, [Chem.MolFromSmiles(r) for r in reactants])
            score, reasons = advanced_evaluation([Chem.MolFromSmiles(r) for r in reactants], products, selected_reaction)

            product_smiles = [Chem.MolToSmiles(prod) for prod in products]
            product_cmls = [smiles_to_cml_with_2d(smiles) for smiles in product_smiles]

            return jsonify({
                "success": True,
                "reactants": reactants,
                "reactions": [{
                    "reaction_name": selected_reaction,
                    "score": score,
                    "products": product_smiles,
                    "products_cmls": product_cmls,
                    "reasons": reasons,
                    "smarts": reaction_smarts
                }],
                "total_reactions_found": 1 if products else 0
            })

        
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
