
import pyreadr
import numpy as np
import pandas as pd
import sys
base_dir = "dir/to/Clipper"
#base_dir = "/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/"
sys.path.append(base_dir)
import Clipper

if __name__ == "__main__":
    exp_d = pyreadr.read_r('/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/data/exp_d.rda')['exp_d'].values
    back_d = pyreadr.read_r('/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/data/back_d.rda')['back_d'].values
    exp_e = pyreadr.read_r('/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/data/exp_e.rda')['exp_e'].values
    back_e = pyreadr.read_r('/Users/gexinzhou/Dropbox/clipper/Github/Clipper_python/data/back_e.rda')['back_e'].values

    print(f"...starting the first test on sample differential analysis")
    re1 = Clipper.clipper(score_exp=exp_d, score_back=back_d, analysis="differential", FDR=[0.01, 0.05, 0.1])
    trueid = np.arange(2000)
    discoveries = pd.DataFrame(re1["discoveries"][0])
    print(f"...the proportion of left-out discoveries: {np.sum(~discoveries.isin(trueid).values)/len(trueid)}")
    print()

    print(f"...starting the second test on sample enrichment analysis")
    re2 = Clipper.clipper(score_exp=exp_e, score_back=back_e, analysis="enrichment", FDR=[0.01, 0.05, 0.1])
    trueid = np.arange(1000)
    discoveries = pd.DataFrame(re2["discoveries"][0])
    print(f"...the accuracy of discoveries: {np.sum(discoveries.isin(trueid).values)/len(trueid)}")
    print()

    print(f"...starting the third test on the sample DEG data")
    DEG_exp_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/DEG_exp_d.RData')['score_exp'].values
    DEG_back_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/DEG_back_d.RData')['score_back'].values
    re_clipper = Clipper.clipper(DEG_exp_d, DEG_back_d, FDR=0.05, analysis="differential")
    print(f"...result is \n{re_clipper}")
    expected_result = """
        {'contrast.score': 'max',
    'contrast.score.value': array([-0.50045148, -2.97313391, -0.92149901, ...,  0.5559539 ,
            -0.8190085 , -1.34792436]),
    'FDR': array([0.05]),
    'contrast.score.thre': array([4.86856861]),
    'discoveries': array([[  675,  2172,  2734,  2837,  3101,  3102,  3104,  3483,  3489,
            3510,  5694,  5704,  7011,  7026,  8957, 12913, 15517, 16257,
            16953, 21321, 22247, 26295, 26893, 28386, 28549, 29490, 29527,
            30817, 31127, 32550, 32691, 32757, 32774, 32833, 32835, 32908,
            34189, 34281, 35583, 36303, 36627, 36636, 38168, 38379, 39083,
            39121, 41044, 41986, 42621, 42622, 42776, 43423, 43674, 44281,
            44704, 44722, 44918, 45027, 46707, 49100, 49129, 51217]]),
    'q': array([1., 1., 1., ..., 1., 1., 1.])}
    """
    print(f"expected results: {expected_result}")
    print()

    print(f"...starting the fourth test on sample Hi_C data")
    Hi_c_exp_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/hi_c_score_exp.rds')[None].values
    Hi_c_back_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/hi_c_score_back.rds')[None].values
    re_hi_c = Clipper.clipper(Hi_c_exp_d, Hi_c_back_d, FDR=[0.01, 0.05, 0.1], analysis="differential")
    Hi_c_exp1 = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/hi_c_exp1.rds')[None].values
    print(f"...Hi_c_exp1[re_hi_c['discoveries'][1]] is \n{Hi_c_exp1[re_hi_c['discoveries'][1]]}")
    expected_result = """
    array([[4.60e+07, 1.11e+08, 4.08e+02],
       [4.70e+07, 1.11e+08, 3.38e+02],
       [4.80e+07, 1.11e+08, 1.98e+02],
       ...,
       [1.97e+08, 2.48e+08, 1.95e+02],
       [1.06e+08, 2.49e+08, 2.30e+01],
       [1.48e+08, 2.49e+08, 1.00e+00]])
    """
    print(f"...expected result: {expected_result}")
    print()

    print(f"...starting the fifth test on sample MACS2 data")
    MACS_exp_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/MACS_exp_d.rds')[None].values
    MACS_back_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/MACS_back_d.rds')[None].values
    re = Clipper.clipper(MACS_exp_d, MACS_back_d, analysis="enrichment")
    macs2_peak = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/macs2_peak.rds')[None]
    def peak_apply(k):
        start, end = macs2_peak.V2[k], macs2_peak.V3[k]
        return np.median(MACS_exp_d[start:end] - MACS_back_d[start:end])
    median_peak = np.array(list(map(peak_apply, np.arange(len(macs2_peak)))))
    clipper_peak = macs2_peak[median_peak >= re["contrast.score.thre"]]
    print(f"...peaks discovered by Clipper: {clipper_peak}")
    print(f"...whether the length is expected: {len(clipper_peak) == 744}")
    print()

    print(f"...starting the sixth test on sample peptide data")
    peptide_exp_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/peptide_score_exp.rds')[None].values
    peptide_back_d = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/peptide_score_back.rds')[None].values
    peptide_match = pyreadr.read_r('/home/wbguo/iproject/Clipper/data/peptide_match.rds')[None].values
    re = Clipper.clipper(peptide_exp_d, peptide_back_d, FDR = 0.01, analysis="enrichment")
    print(f"...peptide_match[re['discoveries'][0]][:6] is \n{peptide_match[re['discoveries'][0]][:6]}")
    expected_result = """
    array([['am190328_006.raw:1572:KYEESGKPR'],
       ['am190328_006.raw:1720:RADSMISDEKER'],
       ['am190328_005.raw:1784:MMVEVAK'],
       ['am190328_006.raw:1807:YPESNYMHK'],
       ['am190328_005.raw:1830:MDMVNYNQK'],
       ['am190328_006.raw:1852:MDMVNYNQK']], dtype=object)
    """
    print(f"...expected result: {expected_result}")
    

