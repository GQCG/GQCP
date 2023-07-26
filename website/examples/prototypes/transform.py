import numpy as np

def tensorToMOBasis(tensor, matrix1, matrix2, matrix3, matrix4):

    temp_1 = np.einsum('as,pqra->pqrs', matrix4, tensor)
    temp_2 = np.einsum('lr,pqls->pqrs', matrix3, temp_1)
    temp_3 = np.einsum('nq,pnrs->pqrs', matrix2, temp_2)
    tensor_nb = np.einsum('mp,mqrs->pqrs', matrix1, temp_3)

    return tensor_nb
