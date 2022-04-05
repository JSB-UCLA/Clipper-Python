import warnings
import numpy as np
from .helper_functions import *

def clipper(score_exp, score_back, analysis, FDR = 0.05, 
            procedure = None, contrast_score = None,
            num_permutation = None, seed = 12345):

    analysis = match_arg(analysis, ["differential", "enrichment"], "analysis")

    score_exp = np.atleast_2d(score_exp)
    score_back = np.atleast_2d(score_back)

    FDR = np.atleast_1d(FDR)
    if (procedure is None):
        procedure = "GZ" if analysis == "differential" else "BC"
    else:
        procedure = match_arg(procedure, ["BC", "aBH", "GZ"], "procedure")

    if (analysis == "differential"):
        contrast_score = match_arg(contrast_score, ["diff", "max"], "contrast_score") if contrast_score is not None else "max"

        re = clipper2sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                        nknockoff=num_permutation, 
                        contrastScore_method= contrast_score, importanceScore_method="diff",
                        FDR_control_method=procedure, ifpowerful=False, seed=seed)
        
        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))
        if np.any(np.array(FDR_nodisc) & (contrast_score == "max")):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold. ")

    elif (analysis == "enrichment"):

        if (np.shape(score_exp)[1] != np.shape(score_back)[1]):
            procedure = "GZ"
        if (contrast_score is None):
            contrast_score = "diff" if procedure == "BC" else "max"
        else:
            contrast_score = match_arg(contrast_score, ["diff", "max"], "contrast_score")
        if procedure == "BC":
            re = clipper1sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                            nknockoff=num_permutation,
                            importanceScore_method=contrast_score,
                            FDR_control_method=procedure, ifpowerful=False, seed=seed)
        if procedure == "GZ":
            re = clipper1sided(score_exp=score_exp, score_back=score_back, FDR=FDR,
                            nknockoff=num_permutation,
                            contrastScore_method=contrast_score,
                            FDR_control_method=procedure, ifpowerful=False, seed=seed)

        FDR_nodisc = np.array(list(map(lambda re_i: len(re_i["discovery"]) == 0, re["results"])))
        if np.any(np.array(FDR_nodisc) & (procedure != 'aBH')):
            warnings.warn(f"At FDR = {', '.join([str(x) for x in FDR[FDR_nodisc]])}, no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold. ")

    contrast_score_value = re["contrastScore"]
    thre = []
    discoveries = []
    results = re["results"]
    for result in results:
        try:
            thre.append(result["thre"])
            discoveries.append(result["discovery"])
        except:
            continue
    thre = np.array(thre)
    discoveries = np.array(discoveries)
    q = results[0]["q"]
    re = {"contrast.score": contrast_score, 
            "contrast.score.value": contrast_score_value,
            "FDR": FDR,
            "contrast.score.thre": thre,
            "discoveries": discoveries,
            "q": q}
    return re