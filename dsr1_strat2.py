from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.AllChem import ReactionFromSmarts
from collections import defaultdict

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

def evaluate_reaction_fit(reactants, products, reaction_name):
    """
    A simple scoring function to evaluate how well a reaction fits.
    Returns a score and reasons.
    """
    if not products:
        return 0.0, ["No products formed"]
    
    score = 0.0
    reasons = []
    
    # Basic score for any product formation
    score += 1.0
    reasons.append("Reaction pattern matched reactants")
    
    # Penalize reactions that generate very small fragments (e.g., water, salts)
    # but keep them if they are expected byproducts
    main_products = []
    for prod in products:
        mol_wt = Descriptors.MolWt(prod)
        if mol_wt > 50:  # Arbitrary threshold
            main_products.append(prod)
            score += 0.5
        else:
            reasons.append(f"Generated small fragment: {Chem.MolToSmiles(prod)}")
    
    if not main_products:
        score -= 2.0
        reasons.append("Only small fragments generated")
    
    # Check if products seem chemically reasonable
    for prod in main_products:
        try:
            # Products should typically have no formal charges or balanced ones
            formal_charge = Chem.GetFormalCharge(prod)
            if formal_charge != 0:
                score -= 0.2
                reasons.append(f"Product has formal charge: {formal_charge}")
        except:
            pass
    
    return max(0.0, score), reasons


def advanced_evaluation(reactants, products, reaction_name):
    """A more sophisticated scoring function."""
    if not products:
        return 0.0, ["No reaction"]
    
    score = 0.0
    reasons = []
    
    # Reaction-specific scoring
    if reaction_name == 'amide_formation':
        # Check if we have the right reactants: acid + amine
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
    
    # Mass should be approximately conserved
    mass_balance_ratio = product_mw / reactant_mw
    if 0.8 < mass_balance_ratio < 1.2:
        score += 1.0
        reasons.append("Good mass balance")
    else:
        score -= 1.0
        reasons.append(f"Poor mass balance: ratio {mass_balance_ratio:.2f}")
    
    # Check for common byproducts
    common_byproducts = {'O', 'CO', 'CCO', 'N', '[NH3+]', '[Cl-]', '[Na+]'}
    product_smiles = {Chem.MolToSmiles(p) for p in products}
    unexpected_small_molecules = product_smiles.intersection(common_byproducts)
    
    if unexpected_small_molecules:
        score -= 0.5
        reasons.append(f"Generated unexpected small molecules: {unexpected_small_molecules}")
    
    return max(0.0, score), reasons




def find_best_reactions(reactant_smiles_list, reaction_library):
    """
    Try all reactions in the library and return scored results.
    """
    reactants = [Chem.MolFromSmiles(smi) for smi in reactant_smiles_list]
    
    results = []
    
    for rxn_name, rxn_smarts in reaction_library.items():
        try:
            products = run_reaction(rxn_smarts, reactants)
            #score, reasons = evaluate_reaction_fit(reactants, products, rxn_name)
            score, reasons = advanced_evaluation(reactants, products, rxn_name)
            
            results.append({
                'reaction_name': rxn_name,
                'score': score,
                'products': products,
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

# Example usage
if __name__ == "__main__":
    # Your reaction library
    reaction_library = {
        'amide_formation': '[C:1](=[O:2])-[OD1:3].[N!H0:4]>>[C:1](=[O:2])[N:4]',
        'esterification': '[C:1](=[O:2])-[OD1:3].[O!H0:4]-[C:5]>>[C:1](=[O:2])[O:4][C:5]',
        'nitration': '[cH:1]>>[c:1][N+](=O)[O-]',
        'hydrogenation': '[C:1]=[C:2]>>[C:1][C:2]',
        'ester_hydrolysis': '[O:1]=[C:2][O:3][C:4]>>[O:1]=[C:2][OH].[OH:3][C:4]',
        'reductive_amination': '[C:1]=[O:2].[N!H0:3]>>[C:1][N:3]',
    }
    
    # Test case: acetic acid + methylamine
    #reactants = ['CC(=O)O', 'CN'] 
    reactants = ['CCC=O', 'CN'] 
    
    results = find_best_reactions(reactants, reaction_library)
    
    print(f"Testing reactants: {' + '.join(reactants)}")
    print("=" * 50)
    
    for result in results:
        if result['score'] > 0:
            print(f"\n{result['reaction_name']} (Score: {result['score']:.2f})")
            print(f"Reasons: {', '.join(result['reasons'])}")
            for i, prod in enumerate(result['products']):
                print(f"  Product {i+1}: {Chem.MolToSmiles(prod)}")




