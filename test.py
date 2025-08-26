from rdkit import Chem
from rdkit.Chem import AllChem

# Define the reactant using SMILES
reactant_smiles = "CCCO"  # Ethanol
reactant = Chem.MolFromSmiles(reactant_smiles)

# Define the reaction using SMARTS
reaction_smarts = "[C:1][OH:2]>>[C:1]=O"  # Oxidation of alcohol to aldehyde
reaction_smarts = "[CH2:1][OH]>>[CH:1]=O"
reaction_smarts = "[CH2,CH3:1]-[OH]>>[C:1]=O"


reaction_smarts = "[cH:1]>>[c:1][N+](=O)[O-]"  # aromatic nitration
reaction_smarts = "[C:1]=[C:2]>>[C:1][C:2]" #alkene hydrogenation
reaction_smarts = "[O:1]=[C:2][O:3][C:4]>>[O:1]=[C:2][OH].[OH:3][C:4]" # ester hydrolysis
reaction_smarts = "[CH2:1][OH]>>[CH:1]=[CH2]" # dehydration of alcohols to alkenes
reaction_smarts = "[C:2]=[O:1]>>[C:2][OH:1]" # reduction of ketone to alcohol

amide_formation = '[C:1](=[O:2])-[OD1:3].[N!H0:4]>>[C:1](=[O:2])[N:4]'
# Explanation:
# [C:1](=[O:2])-[OD1:3] : Carboxylic acid (O with exactly one connection)
# [N!H0:4] : An nitrogen atom with at least one H (i.e., an amine)
# >>[C:1](=[O:2])[N:4] : Forms the amide bond

esterification = '[C:1](=[O:2])-[OD1:3].[O!H0:4]-[C:5]>>[C:1](=[O:2])[O:4][C:5]'
# Explanation:
# [O!H0:4]-[C:5] : An oxygen with at least one H (alcohol) connected to carbon

suzuki_coupling = '[c,Br,I:1]-[B](-O)-O.[c,Br,I:2]-[Cl,Br,I:3]>>[c:1]-[c:2]'
# Note: This is a simplified version. Real SMARTS for Pd-catalyzed reactions
# can be more complex to handle specific halides and boronates.

sn2_reaction = '[C:1]-[Cl,Br,I:2].[O,N,S:3]>>[C:1][O,N,S:3]'
# Explanation:
# Breaks C-Halogen bond, forms C-Nu bond.
reductive_amination = '[C:1]=[O:2].[N!H0:3]>>[C:1][N:3]'
# Note: This simplistically shows the final C-N bond formation, ignoring the imine intermediate and reducing agent.

diels_alder = '[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1-[C:2]=[C:3]-[C:4]1-[C:5]-[C:6]'
# Explanation:
# This creates a 6-membered ring from a diene and a dienophile.
# Mapping atoms is crucial here for the correct connectivity.

epoxidation = '[C:1]=[C:2]>>[C:1]1-[C:2]-O-1'
# Explanation:
# Converts a double bond into a 3-membered epoxide ring.

acid_base = '[O,S,N:1]-[H].[O,N:2]-]>>[O,S,N:1]-.[H]-[O,N:2]'
# Example: Carboxylic acid + Base -> Carboxylate
carboxylate_formation = 'C(=O)[OH].[O-]>>C(=O)[O-].O'


reaction = AllChem.ReactionFromSmarts(reaction_smarts)

# Apply the reaction
products = reaction.RunReactants((reactant,))

# Convert products to SMILES and print them
product_smiles = [Chem.MolToSmiles(product[0]) for product in products]
print("Predicted products:", product_smiles)

