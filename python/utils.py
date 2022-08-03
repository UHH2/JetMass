import uproot
import numpy as np
import awkward as ak


def numpy_to_th2(H,x_edges,y_edges,hist_title="",x_title="",y_title="",add_empty_flow=True):
    
    x_taxis = uproot.writing.identify.to_TAxis(
                                    x_title,x_title,
                                    len(x_edges[:-1]),
                                    x_edges[0],x_edges[-1],
                                    x_edges)
    
    y_taxis = uproot.writing.identify.to_TAxis(
                                    y_title,y_title,
                                    len(y_edges[:-1]),
                                    y_edges[0],y_edges[-1],
                                    y_edges)
    
    def add_2d_padding(A):
        nx,ny = A.shape
        return np.concatenate((np.zeros(ny),A.flatten(),np.zeros(ny))).reshape(nx+2,ny)
    
    if(add_empty_flow):
        H = add_2d_padding(add_2d_padding(H).T).T
    
    th2 = uproot.writing.identify.to_TH2x(hist_title,hist_title,
                                            H.flatten(),
                                            1,1,1,1,
                                            1,1,1,1,
                                            np.array([1.]),
                                            x_taxis,y_taxis)
    return th2


def hist_to_th1(H,hist_name=''):
    x_axis = H.axes[0]
    #take values and variances from hist and 'fill' under- and overflow bins
    values = np.concatenate(([0.],H.values(),[0.])) 
    variances  = np.concatenate(([0.],H.variances(),[0.]))

    
    x_taxis = uproot.writing.identify.to_TAxis(
                                    x_axis.name,x_axis.name,
                                    len(x_axis.edges)-1,
                                    x_axis.edges[0],x_axis.edges[-1],
                                    x_axis.edges)
    
    if hist_name == '':
        hist_name = x_axis.name

    th1 = uproot.writing.identify.to_TH1x(hist_name,hist_name,
                                          values, values.sum(),
                                          values.sum(),variances.sum(),
                                          1.,1.,
                                          variances,
                                          x_taxis
                                          )    
    
    return th1

def np_2d_hist_bin_value(vals,x,y):
    # x_bin = np.digitize([x],vals[1])[0]
    # if(x_bin > vals[1].shape[0]-1):
    #     x_bin = vals[1].shape[0]-1
    # elif(x_bin <= 0):
    #     x_bin = 1

    # y_bin = np.digitize([y],vals[2])[0]
    # if(y_bin > vals[2].shape[0]-1):
    #     y_bin = vals[2].shape[0]-1

    # elif(y_bin <=0):
    #     y_bin = 1
    # return vals[0][x_bin-1,y_bin-1]
    if(isinstance(x,float)):
        x = [x]
    x_bin = get_bin_indx(x,vals[1])
    if(isinstance(y,float)):
        y = [y]
    y_bin = get_bin_indx(y,vals[2])
    print(x_bin,y_bin)
    result = vals[0][x_bin,y_bin]
    if(isinstance(x,ak.Array) and isinstance(y,ak.Array)):
        result = ak.Array(result)
    return result


def get_bin_indx(values,edges):
    indx = np.digitize(ak.Array(values),edges)-1
    indx = ak.where(indx==len(edges)-1,len(edges)-2,indx)
    indx = ak.where(indx<0,0,indx)    
    return indx

def root_th2_bin_value(th2,x,y):
    x_bin = th2.GetXaxis().FindFixBin(x)
    if(x_bin > th2.GetXaxis().GetNbins()):
        x_bin = th2.GetXaxis().GetNbins()
    elif(x_bin <= 0):
        x_bin = 1

    y_bin = th2.GetYaxis().FindFixBin(y)
    if(y_bin > th2.GetYaxis().GetNbins()):
        y_bin = th2.GetYaxis().GetNbins()
    elif(y_bin <=0):
        y_bin = 1

    return th2.GetBinContent(x_bin,y_bin)

def get_2d_hist_bin_content(H,x,y):
    if('TH2' in str(type(H))):
        return root_th2_bin_value(H,x,y)
    elif(isinstance(H,tuple)):
        return np_2d_hist_bin_value(H,x,y)
        



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
