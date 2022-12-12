import ROOT
ROOT.TH1.AddDirectory(False)

class CombineWorkspace:
    def __init__(self, file_path_):
        self.w = self.load_workspace(file_path_)
        self.model_s = self.w.pdf("model_s")
        self.data = self.w.data("data_obs")
        self.data_sets = {d.GetName(): d for d in self.data.split(ROOT.RooCategory("CMS_channel", ""))}
        self.cats = self.data_sets.keys()
        self.obs = self.data_sets[self.cats[0]].get().first()
        self.toys = []
        self.toy_datasets = []

    def Print(self):
        self.data.Print("V")
        print(self.cats)

    def load_toys(self, toys_file):
        # toy_datasets = []
        try:
            _file = ROOT.TFile(toys_file)
            toy_dir = _file.Get("toys")
            for toy in toy_dir.GetListOfKeys():
                if "snapshot" not in toy.GetName():
                    self.toy_datasets.append(toy_dir.Get(toy.GetName()))
        except BaseException as e:
            print(e)
            raise IOError("Not able to load toys from " + toys_file)
        self.toys = [
            {
                t.GetName(): t.binnedClone().createHistogram("toy{}_{}".format(itoy, t.GetName()), self.obs)
                for t in toy.split(ROOT.RooCategory("CMS_channel", ""))
            }
            for itoy, toy in enumerate(self.toy_datasets)
        ]

    def load_workspace(self, workspace_path, workspace_name="w"):
        try:
            workspace_file = ROOT.TFile(workspace_path)
            workspace = workspace_file.Get(workspace_name)
            return workspace
        except BaseException as e:
            print(e)
            raise IOError("Not able to load workspace from " + workspace_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--workspace", type=str, help="path to workspace ROOT-file. Should contain workspace names 'w'."
    )
    parser.add_argument("-t", "--toysFile", type=str, help="path to ROOT-file containing toys.")
    parser.add_argument("-o", "--outfile", type=str, default="toys.root")
    parser.add_argument("-M", "--method", type=str, default="convertToys", choices=["convertToys"])
    args = parser.parse_args()

    ws = CombineWorkspace(args.workspace)

    if args.method == "convertToys":
        ws.load_toys(args.toysFile)
        toys_th1_file = ROOT.TFile(args.outfile, "RECREATE")
        for toy in ws.toys:
            for histname, toy_th1 in toy.items():
                toy_th1.Write()
        toys_th1_file.Close()

