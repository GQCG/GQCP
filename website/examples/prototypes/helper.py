import gqcpy

def getAlphaMatrix(parameters):
    return parameters.expansion().alpha.matrix()

def getBetaMatrix(parameters):
    return parameters.expansion().beta.matrix()

def GQCPStates(state_array):
    par = []
    par_f = []

    for i in range(len(state_array)):
        par.append(gqcpy.UTransformation_d(gqcpy.UTransformationComponent_d(state_array[i].loc[0, "alpha_parameters"]), gqcpy.UTransformationComponent_d(state_array[i].loc[0, "beta_parameters"])))
        par_f.append(gqcpy.UTransformation_d(gqcpy.UTransformationComponent_d(state_array[i].loc[0, "beta_parameters"]), gqcpy.UTransformationComponent_d(state_array[i].loc[0, "alpha_parameters"])))

    return par, par_f
