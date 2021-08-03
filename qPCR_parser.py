import argparse
import numpy as np
import pandas as pd
from collections import defaultdict

def read_input(input, sheetName):
    """
    Parameters: An excel file from qPCR and the sheet name (or number) that which the results are located
    Output: A pandas dataframe with columns "Sample Name", "Target Name", and "CT" and with ouliers filtered

    """
    df = pd.read_excel(input, sheet_name=sheetName)

    # Ensures that the starting row and column of dataframe are the relevent rows and columns.
    for row in range(df.shape[0]):
        for col in range(df.shape[1]):
            if df.iat[row,col] == "Sample Name":
                row_start = row
    df = df.loc[row_start:]
    new_header = df.iloc[0]
    df = df[1:]
    df.columns = new_header
    # Sets the dataframe index to the well number.
    df.set_index("Well", inplace=True)
    df = df.loc[:, df.columns.intersection(['Sample Name','Target Name', "CT", "Ct SD", "Ct Mean"])]
    df = df.dropna()

    # Removes outliers that are greater than the specified number of standard deviations from the mean (default=3)
    df = df[np.abs(df.CT-df.CT.mean()) <= (args.std*df.CT.std())]

    return df

def foldGeneExp(df):
    """
    Parameters: A dataframe that is ready for processing with outliers remvoved.
    Output: A bar chart showing fold gene expressions changes of all samples with repect to control conditions.

    """

    # Turns all samples into list of tuples with format ("sample", "target")
    tuples = [df.filter(["Sample Name", "Target Name"]).to_records(index=False)]
    dictOfMeans = defaultdict(list)

    # list of tuples (Sample, Target)
    sampTargets = []
    samples = []
    targets = []

    # cleans up list of tuples and creates list of samples and list of targets
    for i in tuples:
        for j in range(len(i)):
            sampTargets.append(i[j])
            samples.append(i[j][0])
            targets.append(i[j][1])

    # list of CT values
    CTs = df["CT"].tolist()

    # makes a dictionary of (sample, target) tuples (keys) and their correspoding CT values as a list of values.

    for i in range(len(sampTargets)):
        if i == 0:
            dictOfMeans[samples[0], targets[0]].append(CTs[0])

        elif i != 0 and i != len(sampTargets) - 1:
            # Ensures all samples and targets are included in the dictionary
            if str(sampTargets[i]) == str(sampTargets[i + 1]) and str(sampTargets[i]) == str(sampTargets[i - 1]):
                dictOfMeans[samples[i], targets[i]].append(CTs[i])

            elif str(sampTargets[i]) == str(sampTargets[i + 1]):
                dictOfMeans[samples[i], targets[i]].append(CTs[i])

            elif str(sampTargets[i]) == str(sampTargets[i - 1]) and str(sampTargets[i]) != str(sampTargets[i + 1]):
                dictOfMeans[samples[i], targets[i]].append(CTs[i])
        else:
            dictOfMeans[samples[i], targets[i]].append(CTs[i])


    # Computes the average CT values for sample, target pairs and assigns to a new dictionary
    newDicOfMeans = {}
    dicOfControl = {}

    controls = args.control
    for key, value in dictOfMeans.items():
        if key[0] not in controls:
            newDicOfMeans[key] = sum(value)/len(value)
        else:
            dicOfControl[key] = sum(value)/len(value)

    # Removes duplicates from the list of targets and samples and removes the housekeeper gene from list.

    housekeeper = [args.hg]
    targets1 = list(dict.fromkeys(targets))
    samples1 = list(dict.fromkeys(samples))
    for i in targets1:
        if i in housekeeper:
            targets1.remove(i)

    dictOfHousekeeper_exp = {}
    dictOfTargets_exp = {}
    dictOfTargets_ctrl = {}
    dicOfHousekeeper_ctrl= {}

    for i in housekeeper:
        dictOfHousekeeper_exp[i] = []
    for i in targets1:
        dictOfTargets_exp[i] = []
    for i in targets1:
        dictOfTargets_ctrl[i] = []
    for i in housekeeper:
        dicOfHousekeeper_ctrl[i] = []


    for key, value in newDicOfMeans.items():
        for target in dictOfTargets_exp.keys():
            if str(key[1]) == str(target):
                dictOfTargets_exp[target].append(value)

    for key, value in dicOfControl.items():
        for target in dictOfTargets_ctrl.keys():
            if str(key[1]) == str(target):
                dictOfTargets_ctrl[target].append(value)

    for key, value in newDicOfMeans.items():
        for target in dictOfHousekeeper_exp.keys():
            if str(key[1]) == str(target):
                dictOfHousekeeper_exp[target].append(value)

    for key, value in dicOfControl.items():
        for target in dicOfHousekeeper_ctrl.keys():
            if str(key[1]) == str(target):
                dicOfHousekeeper_ctrl[target].append(value)

        #Creates a new dictionary for each target and computes the delta Ct value for each against the correcsponding housekeeper CT value

    dCT_ctrl = defaultdict(list)

    for key, value in dictOfTargets_ctrl.items():
        for k, v in dicOfHousekeeper_ctrl.items():
            for i in range(len(v)):
                dCT_ctrl[key].append(abs(v[i]-value[i]))


    dCT_exp = defaultdict(list)

    for key, value in dictOfTargets_exp.items():
        for k, v in dictOfHousekeeper_exp.items():
            for i in range(len(v)):
                dCT_exp[key].append(abs(v[i]-value[i]))

    dCT_merged = {key: dCT_ctrl[key] + dCT_exp[key] for key in dCT_ctrl}

    # Computes the average delta CT of all controls - will use this value to compute delta delta CT for all samples
    dCT_ctrl_average = {}
    for key, value in dCT_ctrl.items():
        dCT_ctrl_average[key] = sum(value)/len(value)

    # Creates a new dictionary that compares dCT to the control group
    ddCT = defaultdict(list)

    for key, value in dCT_merged.items():
        for i in range(len(value)):
            ddCT[key].append(value[i] - dCT_ctrl_average[key])

    foldGeneExp = defaultdict(list)

    for key, value in ddCT.items():
        for i in range(len(value)):
            foldGeneExp[key].append(2**(-(value[i])))

    data = pd.DataFrame(foldGeneExp).T
    for column in range(0, len(data)):
        data.rename(columns={column: samples1[column]}, inplace=True)

    ax = data.plot(kind="bar", figsize=(20,10))
    x_offset = -0.03
    y_offset = 0.02
    for p in ax.patches:
        b = p.get_bbox()
        val = "{:+.2f}".format(b.y1 + b.y0)
        ax.annotate(val, ((b.x0 + b.x1)/2 + x_offset, b.y1 + y_offset))
    ax.set_ylabel("Fold gene expression")

    plt.savefig(args.out)

    return None

if __name__ == "__main__":
    p=argparse.ArgumentParser(description="Calculate Fold gene expression and make a bar chart of results from qPCR results excel file")
    p.add_argument("-d", "--data", action="store", dest="data", required=True, help= 'The file you want to analyse')
    p.add_argument("-s", "--sheetName", action="store", dest="sheet", default=0, help= "str, int. the sheet name of your excel file you want to analyse")
    p.add_argument("-hg", '--housekeeper', dest="hg", default=None, help= "Name of your housekeeper gene")
    p.add_argument("-c", "--control", nargs="+", dest="control", required=True, help="the control group name that which you want to compare against")
    p.add_argument("-o", "--output", action="store", default="foo", dest="out", help="output file name to to which you want fold gene expression graph saved to")
    p.add_argument("--std", action="store", type=float,  default=3, dest="std", help="Standard deviation for filtering Ct values.")
    args = p.parse_args()
    foldGeneExp(read_input(args.data, args.sheet))
