from rdkit import Chem
from rdkit.Chem import AllChem
import webbrowser
from rdkit.Chem.Draw import MolToFile
import pubchempy as pcp


print('#############    Organic Final Project   ##############')

print("\nThis simple program prompt the user to input 2 structures (in the SMILES [Simplified Molecular Input Line Entry System] format), outputs the IUPAC names and structures of these reactants, run the reaction between them based on a set of predefined mechanisms, and outputs the final predicted product or products of the reaction in SMILES format, for which the IUPAC name and structure will be also outputted.")

print("\nIt is recommended to place this program file into a foler so that the output files to be generated are stored within it.")

print("\nFive reactions selected from different chapters of the course were included in this program. The selection of these reactions was arbitrary, yet any reaction could be added to the code and the program will be able to predict its products.")

print("\nThese five reactions are:\n\n1. Esterification: input a carboxylic acid as reactant 1 and an alcohol as reactant 2 ; Example: reactant1: 'CC(=O)O' ; reactant2: 'CCO'.")

print("2. Amide preparation: input a carboxylic acid as reactant 1 and an amine as reactant 2 ; Example: reactant1: 'CC(=O)O' ; reactant2: 'CN'.")

print("3. Anhydride preparation: input a carboxylic acid as reactant 1 and an another carboxylic acid as reactant 2 ; Example: reactant1: 'CC(=O)O' ; reactant2: 'CC(=O)O'.")

print("4. Grignard addition to ketone: input a Grignard molecule as reactant 1 and a ketone as reactant 2 ; Example: reactant1: 'C(C)(C)[Mg]Br' ; reactant2: 'CC(=O)C'.")

print("5. Copper reagent addition to aldehyde: input an aldehyde as reactant 1 and a copper reagent as reactant 2 ; Example: reactant1: 'O=CCCC(=O)C' ; reactant2: '[Cu]CC'.")


print("\nIt is recommended to follow the standard convention of placing the more reactive or electrophilic reactant first, and the less reactive or nucleophilic reactant second.")

print("")

print('You might use the reactants examples (given in SMILES format) provided above for each of the 5 reactions for testing purposes.\n')

reactant1 = input('Enter reactant 1 (in SMILES format): ')
reactant2 = input('Enter reactant 2 (in SMILES format): ')
 
input("\nGetting reactants IUPAC names (through PubChem search). Note that if the server is busy, names will be left as unknown. Press Enter to proceed: ")

try:
    name1 = pcp.get_compounds(reactant1, 'smiles',timeout=10)
except TimeoutError:
    print("TimeOut. PubChem server appears to be busy. Try again!")
try:
    name2 = pcp.get_compounds(reactant2, 'smiles',timeout=10)
except TimeoutError:
    print("TimeOut. PubChem server appears to be busy. Try again!")
    
    
n1="unknown1"
n2="unknown2"

if len(name1) > 0:
    name = name1[0].iupac_name
    print(f"The name of the molecule with SMILES '{reactant1}' is '{name}'.")
    n1=name
else:
    print(f"No compound with SMILES '{reactant1}' was found in the PubChem database.")


if len(name2) > 0:
    name = name2[0].iupac_name
    print(f"The name of the molecule with SMILES '{reactant2}' is '{name}'.")
    n2=name
else:
    print(f"No compound with SMILES '{reactant2}' was found in the PubChem database.")


r1 = Chem.MolFromSmiles(reactant1)
r2 = Chem.MolFromSmiles(reactant2)

input("\nTo output reactants structures, press Enter (the structure will be saved as png file in same directory): ")

MolToFile(r1,n1+'.png')
MolToFile(r2,n2+'.png')

reactions=list()


rxn1 = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH].[C:3][OH]>>[C:1](=[O:2])[O][C:3]')
#acid-alcohol (esterification) example : ['CC(=O)O' ; 'CCO']
reactions.append(rxn1)

rxn2 = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]')
#amide-rx examle : ['CC(=O)O' ; 'CN']
reactions.append(rxn2)

rxn3 = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH].[C:3](=[O:4])[OH]>>[C:1](=[O:2])O[C:3](=[O:4])[C]')
#anhydride rx example : ['CC(=O)O' ; 'CC(=O)O'] 
reactions.append(rxn3)

rxn4 = AllChem.ReactionFromSmarts('[Mg][Br].[C:1][C:2](=[O:3])[C:4]>>[C:1][C:2]([MgX:5])[C:4]([OH:3])')
#ketone/grignard rx example : ['C(C)(C)[Mg]Br' ; 'CC(=O)C']
reactions.append(rxn4)

rxn5 = AllChem.ReactionFromSmarts('[C:1]=[O:2].[Cu:3]CC>>[C:1][Cu:3].[O:2]')
#copper-catalyzed conjugate addition example : ['O=CCCC(=O)C' ; '[Cu]CC']
reactions.append(rxn5)


products = reactions[0].RunReactants((r1,r2))
i=0
while len(products) == 0 and i < len(reactions):
    products = reactions[i].RunReactants((r1,r2))
    i=i+1

pro='product'

for product in products[0]:
    print("\nFinal product: "+Chem.MolToSmiles(product))
    try:
        p = pcp.get_compounds(Chem.MolToSmiles(product), 'smiles',timeout=10)
        if not p[0].iupac_name == None:
            pro=p[0].iupac_name
            print(f"The name of the molecule with SMILES '{Chem.MolToSmiles(product)}' is '{pro}'.")
    except TimeoutError:
        print("TimeOut. PubChem server appears to be busy. Product name will be left as product")
    input("\nTo output product structure, press Enter (the structure will be saved as png file in same directory): ")
    MolToFile(product,pro+'.png')
    
input("\nPress Enter to stop execution: ")
