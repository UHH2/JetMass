from utils import bcolors
from coffea import processor


class CoffeaWorkflow(object):
    def __init__(self, name="GenericCoffeaWorkflow"):
        self.workflow_name = name

        import argparse

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--selection", default="vjets")
        self.parser.add_argument("--scaleout", type=int, default=-1)
        self.parser.add_argument("--chunk", type=int, default=-1)
        self.parser.add_argument("--maxchunks", type=int, default=-1)
        self.parser.add_argument("--debug", action="store_true")
        self.parser.add_argument("--megadebug", action="store_true")
        self.args = None

        self.processor_args = {}
        self.processor_instance = None
        self.processor_schema = None

        self.executor = None
        self.executor_args = {}

        self.client = None
        self.cluster = None
        self._init_dask = False

    def parse_args(self):
        self.args = self.parser.parse_args()
        return self.args

    def init_dask_local_client(self):
        if self._init_dask:
            print(
                f"{bcolors.WARNING}WARNING: You did already initialize a dask HTC client!{bcolors.ENDC}"
            )
            return
        from dask.distributed import Client

        self.client = Client()
        self._init_dask = True

    def init_dask_htcondor_client(self, cores=8, memory=24, disk=5):
        if self._init_dask:
            print(
                f"{bcolors.WARNING}WARNING: You did already initialize a dask HTC client!{bcolors.ENDC}"
            )
            return

        from dask.distributed import Client
        from dask_jobqueue import HTCondorCluster

        if self.args.chunk > 0:
            self.processor_args["chunksize"] = self.args.chunk
        if self.args.maxchunks > 0:
            self.processor_args["maxchunks"] = self.args.maxchunks

        if self.processor_instance is None:
            raise (BaseException("You need to set a processor_instance"))

        cluster = HTCondorCluster(cores=cores, memory=f"{memory}GB", disk=f"{disk}GB")
        cluster.scale(self.args.scaleout)
        self.client = Client(cluster)
        self._init_dask = True

    def run_uproot_job_dask(self, samples):

        if self.client is None:
            raise (
                BaseException(
                    "You need to initialize the dask client befor running the job!!"
                )
            )

        self.client.wait_for_workers(1)
        print(f"{self.client}\n{self.client.dashboard_link}")
        # print(self.client.scheduler_info())
        # os.system('ls -lah')
        self.executor = processor.dask_executor
        self.executor_args = {
            "client": self.client,
            "schema": self.processor_schema,
        }
        return self.run_uproot_job(samples)

    def run_uproot_job_local(self, samples):
        self.executor = processor.iterative_executor
        self.executor_args["workers"] = 8
        return self.run_uproot_job(samples)

    def run_uproot_job(self, samples):
        if self.processor_schema is None:
            print(
                bcolors.WARNING,
                "WARNING: You did not specify a processor schema! Falling back to coffea's BaseSchema!!",
                bcolors.ENDC
            )
            from coffea.nanoevents import BaseSchema

            self.procesor_schema = BaseSchema

        output = processor.run_uproot_job(
            samples,
            treename="AnalysisTree",
            processor_instance=self.processor_instance,
            executor=self.executor,
            executor_args=self.executor_args,
            **self.processor_args,
        )
        return output

    def run(self, samples):
        if self.args.debug or self.args.megadebug:
            self.processor_args["maxchunks"] = 1
            self.processor_args["chunksize"] = 1000
            samples = {k: v for k, v in samples.items() if len(v["files"]) > 0}
            for k, v in samples.items():
                samples[k]["files"] = [samples[k]["files"][0]]
            if self.args.megadebug:
                samples = {k: samples[k] for k in list(samples.keys())[2:3]}
                samples["vjets_WJetsMatched"]["files"]=['/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/workdir_vjets_UL18/uhh2.AnalysisModuleRunner.MC.WJetsToQQ_HT800toInf_UL18_10.root']
            print(samples)
        if self.args.scaleout > 0 & self._init_dask:
            return self.run_uproot_job_dask(samples)
        else:
            return self.run_uproot_job_local(samples)
