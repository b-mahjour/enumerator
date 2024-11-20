from rdkit import Chem
from rdkit.Chem import AllChem


def return33():
    A_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6][C:5])=O)") #A - B
    A_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6][C:5])O)") #A - B
    A_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6][C:5]))") #A - B

    A_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][N:10][C:6][C:5])=O)") #A - B
    A_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][N:10][C:6][C:5])O)") #A - B
    A_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([N:10][C:6][C:5]))") #A - B
  
    A_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([N:10][C:6]([C:5]))") #alpha - B

    A_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][N:10][C:6]([C:5])") #beta - B


    alpha_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5])[N:10])=O)") #A - alpha
    alpha_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5])[N:10])O)") #A - alpha
    alpha_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5])[N:10]))") #A - alpha

    alpha_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6]([C:5])[N:10])=O)") #A - alpha
    alpha_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6]([C:5])[N:10])O)") #A - alpha
    alpha_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6]([C:5])[N:10]))") #A - alpha

    alpha_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:6]([C:5])([N:10]))") #alpha - alpha

    alpha_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:6]([C:5])([N:10])") #beta - alpha
   

    beta_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6]([N:10]))=O)") #A - beta
    beta_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6]([N:10]))O)") #A - beta
    beta_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6]([N:10])))") #A - beta

    beta_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5][C:6]([N:10]))=O)") #A - beta
    beta_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5][C:6]([N:10]))O)") #A - beta
    beta_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5][C:6]([N:10])))") #A - beta

    beta_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:5][C:6][N:10])") #alpha - beta

    beta_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:5][C:6][N:10]") #beta - beta


    alpha_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5]))=O)") #A - alpha
    alpha_mA_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5]))O)") #A - alpha
    alpha_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6]([C:5])))") #A - alpha

    alpha_mA_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6]([C:5]))=O)") #A - alpha
    alpha_mA_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6]([C:5]))O)") #A - alpha
    alpha_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([C:6]([C:5])))") #A - alpha
    
    alpha_mA_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:6]([C:5]))") #alpha - alpha
    
    alpha_mA_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:6]([C:5])") #beta - alpha


    beta_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6])=O)") #A - beta
    beta_mA_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6])O)") #A - beta
    beta_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5][C:6]))") #A - beta

    beta_mA_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5][C:6])=O)") #A - beta
    beta_mA_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5][C:6])O)") #A - beta
    beta_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([C:5][C:6]))") #A - beta

    beta_mA_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:5][C:6])") #alpha - beta
    beta_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:5][C:6])") #alpha - beta
    beta_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:5][C:6])") #alpha - beta
    beta_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:5][C:6])") #alpha - beta
    beta_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:5][C:6])") #alpha - beta

    beta_mA_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:5][C:6]") #beta - beta

    all_rxns = [
      A_BC_mOH, A_BC_H, A_BC_HH, A_BO, A_BO_H, A_BO_HH, A_alpha, A_alpha_mOH, A_alpha_H, A_alpha_HH, A_alpha_mB, A_beta, A_beta_mOH, A_beta_H, A_beta_HH, A_beta_mB,
      alpha_BC_mOH, alpha_BC_H, alpha_BC_HH, alpha_BO, alpha_BO_H, alpha_BO_HH, alpha_alpha, alpha_alpha_mOH, alpha_alpha_H, alpha_alpha_HH, alpha_alpha_mB, alpha_beta, alpha_beta_mOH, alpha_beta_H, alpha_beta_HH, alpha_beta_mB,
      beta_BC_mOH, beta_BC_H, beta_BC_HH, beta_BO, beta_BO_H, beta_BO_HH, beta_alpha, beta_alpha_mOH, beta_alpha_H, beta_alpha_HH, beta_alpha_mB, beta_beta, beta_beta_mOH, beta_beta_H, beta_beta_HH, beta_beta_mB,
      alpha_mA_BC_mOH, alpha_mA_BC_H, alpha_mA_BC_HH, alpha_mA_BO, alpha_mA_BO_H, alpha_mA_BO_HH, alpha_mA_alpha, alpha_mA_alpha_mOH, alpha_mA_alpha_H, alpha_mA_alpha_HH, alpha_mA_alpha_mB, alpha_mA_beta, alpha_mA_beta_mOH, alpha_mA_beta_H, alpha_mA_beta_HH, alpha_mA_beta_mB,
      beta_mA_BC_mOH, beta_mA_BC_H, beta_mA_BC_HH, beta_mA_BO, beta_mA_BO_H, beta_mA_BO_HH, beta_mA_alpha, beta_mA_alpha_mOH, beta_mA_alpha_H, beta_mA_alpha_HH, beta_mA_alpha_mB, beta_mA_beta, beta_mA_beta_mOH, beta_mA_beta_H, beta_mA_beta_HH, beta_mA_beta_mB
    ]

    all_rxns_names = [
      "A_BC_mOH", "A_BC_H", "A_BC_HH", "A_BO", "A_BO_H", "A_BO_HH", "A_alpha", "A_alpha_mOH", "A_alpha_H", "A_alpha_HH", "A_alpha_mB", "A_beta", "A_beta_mOH", "A_beta_H", "A_beta_HH", "A_beta_mB",
      "alpha_BC_mOH", "alpha_BC_H", "alpha_BC_HH", "alpha_BO", "alpha_BO_H", "alpha_BO_HH", "alpha_alpha", "alpha_alpha_mOH", "alpha_alpha_H", "alpha_alpha_HH", "alpha_alpha_mB", "alpha_beta", "alpha_beta_mOH", "alpha_beta_H", "alpha_beta_HH", "alpha_beta_mB",
      "beta_BC_mOH", "beta_BC_H", "beta_BC_HH", "beta_BO", "beta_BO_H", "beta_BO_HH", "beta_alpha", "beta_alpha_mOH", "beta_alpha_H", "beta_alpha_HH", "beta_alpha_mB", "beta_beta", "beta_beta_mOH", "beta_beta_H", "beta_beta_HH", "beta_beta_mB",
      "alpha_mA_BC_mOH", "alpha_mA_BC_H", "alpha_mA_BC_HH", "alpha_mA_BO", "alpha_mA_BO_H", "alpha_mA_BO_HH", "alpha_mA_alpha", "alpha_mA_alpha_mOH", "alpha_mA_alpha_H", "alpha_mA_alpha_HH", "alpha_mA_alpha_mB", "alpha_mA_beta", "alpha_mA_beta_mOH", "alpha_mA_beta_H", "alpha_mA_beta_HH", "alpha_mA_beta_mB",
      "beta_mA_BC_mOH", "beta_mA_BC_H", "beta_mA_BC_HH", "beta_mA_BO", "beta_mA_BO_H", "beta_mA_BO_HH", "beta_mA_alpha", "beta_mA_alpha_mOH", "beta_mA_alpha_H", "beta_mA_alpha_HH", "beta_mA_alpha_mB", "beta_mA_beta", "beta_mA_beta_mOH", "beta_mA_beta_H", "beta_mA_beta_HH", "beta_mA_beta_mB"
    ]

    return all_rxns

def return23():
    A_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6]=[C:5])=O)") #A - B
    A_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6]=[C:5])O)") #A - B
    A_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([N:10][C:6]=[C:5]))") #A - B

    A_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][N:10][C:6]=[C:5])=O)") #A - B
    A_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][N:10][C:6]=[C:5])O)") #A - B
    A_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([N:10][C:6]=[C:5]))") #A - B
  
    A_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([N:10][C:6](=[C:5]))") #alpha - B

    A_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][N:10][C:6](=[C:5])") #beta - B


    alpha_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5])[N:10])=O)") #A - alpha
    alpha_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5])[N:10])O)") #A - alpha
    alpha_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5])[N:10]))") #A - alpha

    alpha_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6](=[C:5])[N:10])=O)") #A - alpha
    alpha_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6](=[C:5])[N:10])O)") #A - alpha
    alpha_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6](=[C:5])[N:10]))") #A - alpha

    alpha_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:6](=[C:5])([N:10]))") #alpha - alpha

    alpha_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:6](=[C:5])([N:10])") #beta - alpha
   

    beta_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6]([N:10]))=O)") #A - beta
    beta_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6]([N:10]))O)") #A - beta
    beta_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6]([N:10])))") #A - beta

    beta_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5]=[C:6]([N:10]))=O)") #A - beta
    beta_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5]=[C:6]([N:10]))O)") #A - beta
    beta_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5]=[C:6]([N:10])))") #A - beta

    beta_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:5]=[C:6][N:10])") #alpha - beta

    beta_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:5]=[C:6][N:10]") #beta - beta


    alpha_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5]))=O)") #A - alpha
    alpha_mA_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5]))O)") #A - alpha
    alpha_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:6](=[C:5])))") #A - alpha

    alpha_mA_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6](=[C:5]))=O)") #A - alpha
    alpha_mA_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:6](=[C:5]))O)") #A - alpha
    alpha_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([C:6](=[C:5])))") #A - alpha
    
    alpha_mA_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:6](=[C:5]))") #alpha - alpha
    
    alpha_mA_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:6](=[C:5])") #beta - alpha


    beta_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6])=O)") #A - beta
    beta_mA_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6])O)") #A - beta
    beta_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([C:5]=[C:6]))") #A - beta

    beta_mA_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5]=[C:6])=O)") #A - beta
    beta_mA_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4][C:5]=[C:6])O)") #A - beta
    beta_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]O([C:5]=[C:6]))") #A - beta

    beta_mA_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4])=O)([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]=O)([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8]([O:4]))([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:8])([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:1][C:2]([C:5]=[C:6])") #alpha - beta

    beta_mA_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2]([C:8])[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1][C:2]([C:8]([O:4])=O)>>[C:2][C:1][C:5]=[C:6]") #beta - beta

    all_rxns = [
      A_BC_mOH, A_BC_H, A_BC_HH, A_BO, A_BO_H, A_BO_HH, A_alpha, A_alpha_mOH, A_alpha_H, A_alpha_HH, A_alpha_mB, A_beta, A_beta_mOH, A_beta_H, A_beta_HH, A_beta_mB,
      alpha_BC_mOH, alpha_BC_H, alpha_BC_HH, alpha_BO, alpha_BO_H, alpha_BO_HH, alpha_alpha, alpha_alpha_mOH, alpha_alpha_H, alpha_alpha_HH, alpha_alpha_mB, alpha_beta, alpha_beta_mOH, alpha_beta_H, alpha_beta_HH, alpha_beta_mB,
      beta_BC_mOH, beta_BC_H, beta_BC_HH, beta_BO, beta_BO_H, beta_BO_HH, beta_alpha, beta_alpha_mOH, beta_alpha_H, beta_alpha_HH, beta_alpha_mB, beta_beta, beta_beta_mOH, beta_beta_H, beta_beta_HH, beta_beta_mB,
      alpha_mA_BC_mOH, alpha_mA_BC_H, alpha_mA_BC_HH, alpha_mA_BO, alpha_mA_BO_H, alpha_mA_BO_HH, alpha_mA_alpha, alpha_mA_alpha_mOH, alpha_mA_alpha_H, alpha_mA_alpha_HH, alpha_mA_alpha_mB, alpha_mA_beta, alpha_mA_beta_mOH, alpha_mA_beta_H, alpha_mA_beta_HH, alpha_mA_beta_mB,
      beta_mA_BC_mOH, beta_mA_BC_H, beta_mA_BC_HH, beta_mA_BO, beta_mA_BO_H, beta_mA_BO_HH, beta_mA_alpha, beta_mA_alpha_mOH, beta_mA_alpha_H, beta_mA_alpha_HH, beta_mA_alpha_mB, beta_mA_beta, beta_mA_beta_mOH, beta_mA_beta_H, beta_mA_beta_HH, beta_mA_beta_mB
    ]
    return all_rxns

def return22():
    A_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6]=[C:5])=O)") #A - B
    A_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6]=[C:5])O)") #A - B
    A_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6]=[C:5]))") #A - B

    A_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][N:10][C:6]=[C:5])=O)") #A - B
    A_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][N:10][C:6]=[C:5])O)") #A - B
    A_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([N:10][C:6]=[C:5]))") #A - B
  
    A_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([N:10][C:6](=[C:5]))") #alpha - B
    A_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([N:10][C:6](=[C:5]))") #alpha - B

    A_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][N:10][C:6](=[C:5])") #beta - B
    A_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][N:10][C:6](=[C:5])") #beta - B


    alpha_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5])[N:10])=O)") #A - alpha
    alpha_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5])[N:10])O)") #A - alpha
    alpha_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5])[N:10]))") #A - alpha

    alpha_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6](=[C:5])[N:10])=O)") #A - alpha
    alpha_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6](=[C:5])[N:10])O)") #A - alpha
    alpha_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6](=[C:5])[N:10]))") #A - alpha

    alpha_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:6](=[C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:6](=[C:5])([N:10]))") #alpha - alpha

    alpha_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:6](=[C:5])([N:10])") #beta - alpha
    alpha_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:6](=[C:5])([N:10])") #beta - alpha
   

    beta_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6]([N:10]))=O)") #A - beta
    beta_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6]([N:10]))O)") #A - beta
    beta_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6]([N:10])))") #A - beta

    beta_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5]=[C:6]([N:10]))=O)") #A - beta
    beta_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5]=[C:6]([N:10]))O)") #A - beta
    beta_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5]=[C:6]([N:10])))") #A - beta

    beta_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:5]=[C:6][N:10])") #alpha - beta
    beta_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:5]=[C:6][N:10])") #alpha - beta

    beta_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:5]=[C:6][N:10]") #beta - beta
    beta_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:5]=[C:6][N:10]") #beta - beta


    alpha_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5]))=O)") #A - alpha
    alpha_mA_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5]))O)") #A - alpha
    alpha_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6](=[C:5])))") #A - alpha

    alpha_mA_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6](=[C:5]))=O)") #A - alpha
    alpha_mA_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6](=[C:5]))O)") #A - alpha
    alpha_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([C:6](=[C:5])))") #A - alpha
    
    alpha_mA_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:6](=[C:5]))") #alpha - alpha
    alpha_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:6](=[C:5]))") #alpha - alpha
    
    alpha_mA_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:6](=[C:5])") #beta - alpha    
    alpha_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:6](=[C:5])") #beta - alpha


    beta_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6])=O)") #A - beta
    beta_mA_BC_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6])O)") #A - beta
    beta_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5]=[C:6]))") #A - beta

    beta_mA_BO = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5]=[C:6])=O)") #A - beta
    beta_mA_BO_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5]=[C:6])O)") #A - beta
    beta_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([C:5]=[C:6]))") #A - beta

    beta_mA_alpha = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:5]=[C:6])") #alpha - beta
    beta_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:5]=[C:6])") #alpha - beta

    beta_mA_beta = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_H = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:5]=[C:6]") #beta - beta
    beta_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5]=[C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:5]=[C:6]") #beta - beta

    all_rxns = [
      A_BC_mOH, A_BC_H, A_BC_HH, A_BO, A_BO_H, A_BO_HH, A_alpha, A_alpha_mOH, A_alpha_H, A_alpha_HH, A_alpha_mB, A_beta, A_beta_mOH, A_beta_H, A_beta_HH, A_beta_mB,
      alpha_BC_mOH, alpha_BC_H, alpha_BC_HH, alpha_BO, alpha_BO_H, alpha_BO_HH, alpha_alpha, alpha_alpha_mOH, alpha_alpha_H, alpha_alpha_HH, alpha_alpha_mB, alpha_beta, alpha_beta_mOH, alpha_beta_H, alpha_beta_HH, alpha_beta_mB,
      beta_BC_mOH, beta_BC_H, beta_BC_HH, beta_BO, beta_BO_H, beta_BO_HH, beta_alpha, beta_alpha_mOH, beta_alpha_H, beta_alpha_HH, beta_alpha_mB, beta_beta, beta_beta_mOH, beta_beta_H, beta_beta_HH, beta_beta_mB,
      alpha_mA_BC_mOH, alpha_mA_BC_H, alpha_mA_BC_HH, alpha_mA_BO, alpha_mA_BO_H, alpha_mA_BO_HH, alpha_mA_alpha, alpha_mA_alpha_mOH, alpha_mA_alpha_H, alpha_mA_alpha_HH, alpha_mA_alpha_mB, alpha_mA_beta, alpha_mA_beta_mOH, alpha_mA_beta_H, alpha_mA_beta_HH, alpha_mA_beta_mB,
      beta_mA_BC_mOH, beta_mA_BC_H, beta_mA_BC_HH, beta_mA_BO, beta_mA_BO_H, beta_mA_BO_HH, beta_mA_alpha, beta_mA_alpha_mOH, beta_mA_alpha_H, beta_mA_alpha_HH, beta_mA_alpha_mB, beta_mA_beta, beta_mA_beta_mOH, beta_mA_beta_H, beta_mA_beta_HH, beta_mA_beta_mB
    ]
    return all_rxns

def return32():
    A_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6][C:5])=O)") #A - B
    A_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6][C:5])O)") #A - B
    A_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([N:10][C:6][C:5]))") #A - B

    A_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][N:10][C:6][C:5])=O)") #A - B
    A_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][N:10][C:6][C:5])O)") #A - B
    A_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([N:10][C:6][C:5]))") #A - B
  
    A_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([N:10][C:6]([C:5]))") #alpha - B
    A_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([N:10][C:6]([C:5]))") #alpha - B

    A_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][N:10][C:6]([C:5])") #beta - B
    A_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][N:10][C:6]([C:5])") #beta - B


    alpha_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5])[N:10])=O)") #A - alpha
    alpha_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5])[N:10])O)") #A - alpha
    alpha_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5])[N:10]))") #A - alpha

    alpha_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6]([C:5])[N:10])=O)") #A - alpha
    alpha_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6]([C:5])[N:10])O)") #A - alpha
    alpha_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6]([C:5])[N:10]))") #A - alpha

    alpha_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:6]([C:5])([N:10]))") #alpha - alpha
    alpha_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:6]([C:5])([N:10]))") #alpha - alpha

    alpha_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:6]([C:5])([N:10])") #beta - alpha
    alpha_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:6]([C:5])([N:10])") #beta - alpha
   

    beta_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6]([N:10]))=O)") #A - beta
    beta_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6]([N:10]))O)") #A - beta
    beta_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6]([N:10])))") #A - beta

    beta_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5][C:6]([N:10]))=O)") #A - beta
    beta_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5][C:6]([N:10]))O)") #A - beta
    beta_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5][C:6]([N:10])))") #A - beta

    beta_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:5][C:6][N:10])") #alpha - beta
    beta_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:5][C:6][N:10])") #alpha - beta

    beta_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:5][C:6][N:10]") #beta - beta
    beta_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:5][C:6][N:10]") #beta - beta


    alpha_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5]))=O)") #A - alpha
    alpha_mA_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5]))O)") #A - alpha
    alpha_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:6]([C:5])))") #A - alpha

    alpha_mA_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6]([C:5]))=O)") #A - alpha
    alpha_mA_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:6]([C:5]))O)") #A - alpha
    alpha_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([C:6]([C:5])))") #A - alpha
    
    alpha_mA_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:6]([C:5]))") #alpha - alpha
    alpha_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:6]([C:5]))") #alpha - alpha
    
    alpha_mA_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:6]([C:5])") #beta - alpha    
    alpha_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:6]([C:5])") #beta - alpha


    beta_mA_BC_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6])=O)") #A - beta
    beta_mA_BC_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6])O)") #A - beta
    beta_mA_BC_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([C:5][C:6]))") #A - beta

    beta_mA_BO = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5][C:6])=O)") #A - beta
    beta_mA_BO_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4][C:5][C:6])O)") #A - beta
    beta_mA_BO_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]O([C:5][C:6]))") #A - beta

    beta_mA_alpha = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4])=O)([C:5][C:6])") #alpha - beta
    beta_mA_alpha_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]=O)([C:5][C:6])") #alpha - beta
    beta_mA_alpha_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8]([O:4]))([C:5][C:6])") #alpha - beta
    beta_mA_alpha_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:8])([C:5][C:6])") #alpha - beta
    beta_mA_alpha_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:1]=[C:2]([C:5][C:6])") #alpha - beta

    beta_mA_beta = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4])=O)=[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_mOH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]=O)=[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_H = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8]([O:4]))=[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_HH = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]([C:8])=[C:1][C:5][C:6]") #beta - beta
    beta_mA_beta_mB = AllChem.ReactionFromSmarts("[C:5][C:6]([N:10]).[C:1]=[C:2]([C:8]([O:4])=O)>>[C:2]=[C:1][C:5][C:6]") #beta - beta

    all_rxns = [
      A_BC_mOH, A_BC_H, A_BC_HH, A_BO, A_BO_H, A_BO_HH, A_alpha, A_alpha_mOH, A_alpha_H, A_alpha_HH, A_alpha_mB, A_beta, A_beta_mOH, A_beta_H, A_beta_HH, A_beta_mB,
      alpha_BC_mOH, alpha_BC_H, alpha_BC_HH, alpha_BO, alpha_BO_H, alpha_BO_HH, alpha_alpha, alpha_alpha_mOH, alpha_alpha_H, alpha_alpha_HH, alpha_alpha_mB, alpha_beta, alpha_beta_mOH, alpha_beta_H, alpha_beta_HH, alpha_beta_mB,
      beta_BC_mOH, beta_BC_H, beta_BC_HH, beta_BO, beta_BO_H, beta_BO_HH, beta_alpha, beta_alpha_mOH, beta_alpha_H, beta_alpha_HH, beta_alpha_mB, beta_beta, beta_beta_mOH, beta_beta_H, beta_beta_HH, beta_beta_mB,
      alpha_mA_BC_mOH, alpha_mA_BC_H, alpha_mA_BC_HH, alpha_mA_BO, alpha_mA_BO_H, alpha_mA_BO_HH, alpha_mA_alpha, alpha_mA_alpha_mOH, alpha_mA_alpha_H, alpha_mA_alpha_HH, alpha_mA_alpha_mB, alpha_mA_beta, alpha_mA_beta_mOH, alpha_mA_beta_H, alpha_mA_beta_HH, alpha_mA_beta_mB,
      beta_mA_BC_mOH, beta_mA_BC_H, beta_mA_BC_HH, beta_mA_BO, beta_mA_BO_H, beta_mA_BO_HH, beta_mA_alpha, beta_mA_alpha_mOH, beta_mA_alpha_H, beta_mA_alpha_HH, beta_mA_alpha_mB, beta_mA_beta, beta_mA_beta_mOH, beta_mA_beta_H, beta_mA_beta_HH, beta_mA_beta_mB
    ]
    return all_rxns
