from __future__ import print_function
import os


class GenericJob(object):
    def __init__(self, name, workdir):
        self._workdir = workdir
        self._name = name
        self._wrapper_lines = []
        self._arguments = []
        self._queue = 1
        self.runtime = 10  # hours
        self.runtime *= 60 * 60  # seconds
        self.memory = 4
        self.disk = 4

    def __str__(self):
        return "name: {}, workdir: {}".format(self._name, self._workdir)

    def submit_job(self, debug=False):
        if len(self._wrapper_lines) == 0:
            raise RuntimeError("job wrapper empty!")
        self.prep_job()
        debug_str = " --dry-run {}/{}_dryrun.log ".format(self._workdir, self._name) if debug else ""
        os.system("condor_submit {}/{}.submit{}".format(self._workdir, self._name, debug_str))

    def prep_job(self):
        if not os.path.exists(self._workdir):
            os.makedirs(self._workdir)
        with open("{}/{}_wrapper.sh".format(self._workdir, self._name), "w") as wrapper_file:
            for line in self._wrapper_lines:
                wrapper_file.write("{}\n".format(line))
        os.system("chmod +x {}/{}_wrapper.sh".format(self._workdir, self._name))
        with open("{}/{}.submit".format(self._workdir, self._name), "w") as submit_file:
            submit_file.write('Requirements      = ( OpSysAndVer == "CentOS7" )\n')
            submit_file.write("universe          = vanilla\n")
            submit_file.write("notification      = Error\n")
            submit_file.write("initialdir        = {}\n".format(self._workdir))
            submit_file.write("output            = {}.o$(ClusterId).$(Process)\n".format(self._name))
            submit_file.write("error             = {}.e$(ClusterId).$(Process)\n".format(self._name))
            submit_file.write("log               = {}.$(Cluster).log\n".format(self._name))
            submit_file.write("+RequestRuntime   = {}\n".format(self.runtime))
            submit_file.write("RequestMemory     = {}G\n".format(self.memory))
            submit_file.write("RequestDisk       = {}G\n".format(self.disk))
            submit_file.write("JobBatchName      = {}\n".format(self._name))
            submit_file.write("executable        = {}/{}_wrapper.sh\n".format(self._workdir, self._name))
            if len(self._arguments) > 0:
                submit_file.write("arguments         = {}\n".format(" ".join(self._arguments)))
            submit_file.write("queue {}\n".format(self._queue))


class JobSubmitter(object):
    def __init__(self, base_dir="./"):
        self._base_dir = base_dir
        self._jobs = []

    def add_job(self, name, subworkdir, wrapper_lines):
        self._jobs.append(GenericJob(name, "{}/{}".format(self._base_dir, subworkdir)))
        self._jobs[-1]._wrapper_lines = wrapper_lines

    def add_job_array(self, name, subworkdir, wrapper_lines, arguments, queue):
        self._jobs.append(GenericJob(name, "{}/{}".format(self._base_dir, subworkdir)))
        self._jobs[-1]._wrapper_lines = wrapper_lines
        self._jobs[-1]._arguments = arguments
        self._jobs[-1]._queue = queue

    def submit_jobs(self, debug=False):
        print("submitting {} jobs.".format(len(self._jobs)))
        for job in self._jobs:
            print(job)
            job.submit_job(debug=debug)
