import ROOT

class CombineWorkspace:
    def __init__(self, file_path_):
        self.w = self.load_workspace(file_path_)
        self.model_s = self.w.pdf("model_s")
        self.data = self.w.data("data_obs")
        self.data_sets = {d.GetName():d for d in self.data.split(ROOT.RooCategory("CMS_channel", ""))}
        self.cats = self.data_sets.keys()
        self.obs = self.data_sets[self.cats[0]].get().first()
        self.toys = []
        
    def Print(self):
        self.data.Print("V")
        print(self.cats)

    def load_toys(self, toys_file):
        toys = []
        try:
            _file = ROOT.TFile(toys_file)
            toy_dir = _file.Get("toys")
            for l in toy_dir.GetListOfKeys():
                if "snapshot" not in l.GetName():
                    toys.append(toy_dir.Get(l.GetName()))
        except:
            raise IOError("Not able to load toys from " + toys_file)
        self.toys_sets = toys
        self.toys = [{t.GetName():t for t in toy.split(ROOT.RooCategory("CMS_channel",""))} for toy in toys]
                
    def load_workspace(self, workspace_path,workspace_name='w'):
        workspace = None
        try:
            _file = ROOT.TFile(workspace_path)
            workspace = _file.Get(workspace_name)
        except:
            raise IOError("Not able to load workspace from " + workspace_path)
    
        return workspace
