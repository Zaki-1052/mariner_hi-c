# ML/cmpts/repos/SNIPER/utilities/mm10_config.py
# mm10 (GRCm38) chromosome constants for SNIPER adaptation.

ODD_CHROMS  = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
EVEN_CHROMS = [2, 4, 6, 8, 10, 12, 14, 16, 18]

N_ODD_BINS  = 12808
N_EVEN_BINS = 11831

RESOLUTION  = 100000
SIZES_FILE  = 'data/mm10.chrom.sizes'

N_CLASSES   = 4
SUBC_NAMES  = ['A1', 'A2', 'B1', 'B2']
SUBC_COLORS = ['34,139,34', '152,251,152', '220,20,60', '255,255,0']

TRAIN_VAL_SPLIT = 7000

JUICER_CONTAINER = '/cm/shared/apps/containers/singularity/juicer/juicer_2.0.1.sif'
JUICER_JAR       = '/opt/juicer/CPU/common/juicer_tools.jar'

ODD_BIN_COUNTS  = [1955, 1601, 1519, 1455, 1246, 1221, 1205, 1041, 950, 615]
EVEN_BIN_COUNTS = [1822, 1566, 1498, 1295, 1307, 1202, 1250, 983, 908]
ODD_OFFSETS     = [0, 1955, 3556, 5075, 6530, 7776, 8997, 10202, 11243, 12193]
EVEN_OFFSETS    = [0, 1822, 3388, 4886, 6181, 7488, 8690, 9940, 10923]
