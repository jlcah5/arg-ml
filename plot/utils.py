import pickle
import numpy as np
def getf1(precision, recall):
    """
    Description: given an array with precision and recall at various thresholds, return best f1 score
    Arguments:
        precision: array of precision values
        recall: array of recall values
    """
    for idx, j in enumerate(range(len(precision))):
        best = 0
        bestIdx=0
        if precision[j]+recall[j]==0:
            f1=0
        else:
            f1 = 2*precision[j]*recall[j]/(precision[j]+recall[j])
        if f1 > best:
            best = f1
            bestIdx = idx
        return best


def getPR_size(model, model_type, growth_type, train_size, start=100):
    """
    Description: return average f1, standard error f1, and precision and recall list of a given model
    Arguments:
        model: demographic model
        model_type: neural network type (transformer, cnn)
        growth_type: (growth, noGrowth)
        train_size: number of training samples in kb, in increments of 50k
        start: index to load data
    """
    f1_list=[]
    precision_list = []
    recall_list=[]
    for i in range(start,start+10):
        with open(f'pr_pickle/200kb_0_precision_recall_{model}_{model_type}_{growth_type}_{i}_{train_size}K.pickle', "rb") as f:
            coord = pickle.load(f)
            precision_list.append(np.concatenate([coord[0], [0]*(21-len(coord[0]))]))
            recall_list.append(np.concatenate([coord[1], [0]*(21-len(coord[1]))]))
            best = 0
            bestIidx=0
            for idx, j in enumerate(range(len(coord[0]))):
                if coord[0][j]+coord[1][j]==0:
                    f1=0
                else:
                    f1 = 2*coord[0][j]*coord[1][j]/(coord[0][j]+coord[1][j])
                if f1 > best:
                    best = f1
                    bestIdx = idx
        f1_list.append(best)
    f1_list = np.array(f1_list)
    precision_list = np.vstack(precision_list).sum(axis=0)/10
    recall_list= np.vstack(recall_list).sum(axis=0)/10
    return f1_list.mean(), f1_list.std()/np.sqrt(10), precision_list, recall_list


def getPR_sstar(growth_type, eval_type):
    """
    Description: return average f1, standard error f1, and precision and recall lsit of given sstar experiment
    Arguments:
        growth_type: (growth, noGrowth)
        eval_type: (in, out)
    """
    f1_list=[]
    precision_list = []
    recall_list=[]
    if eval_type == "in":
        start = 4
    else:
        start = 100
    for i in range(start,start+10):
        with open(f'pr_pickle/200kb_0_precision_recall_sstar_{growth_type}_{i}.pickle', "rb") as f:
            coord = pickle.load(f)
            precision_list.append(np.concatenate([coord[0], [0]*(21-len(coord[0]))]))
            recall_list.append(np.concatenate([coord[1], [0]*(21-len(coord[1]))]))
            best = 0
            bestIidx=0
            for idx, j in enumerate(range(len(coord[0]))):
                if coord[0][j]+coord[1][j]==0:
                    f1=0
                else:
                    f1 = 2*coord[0][j]*coord[1][j]/(coord[0][j]+coord[1][j])
                if f1 > best:
                    best = f1
                    bestIdx = idx
        f1_list.append(best)
    f1_list = np.array(f1_list)
    precision_list = np.vstack(precision_list).sum(axis=0)/10
    recall_list= np.vstack(recall_list).sum(axis=0)/10
    return f1_list.mean(), f1_list.std()/np.sqrt(10), precision_list, recall_list