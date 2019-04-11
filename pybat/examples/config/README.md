# Configuring pybat for workflows via fireworks

In this directory you can find some example configuration files for the pybat workflow setup. These files correspond to the typical [Fireworks](https://materialsproject.github.io/fireworks/index.html) YAML files that describe both the fireworker and the launchpad, as well as some examples of jobscript templates for running Fireworks on a queueing system. Here's a brief overview:

1. `launchpad/example_launchpad.yaml` is an example of the launchpad file, i.e. the file that contains all the information about your mongoDB database which is used to store the workflows. Once you've filled in the details of your mongoDB database, you can use `pybat config launchpad -l <launchpad_file>` to configure the launchpad.

2. `fworker/example_fworker.yaml` provides information on the fireworker, i.e. the machine that performs the calculations. Important here is that pybat uses the fireworker `category` to make sure that Fireworks that require a specific amount of nodes are only picked up by a job running on that number of nodes. Also, the `pybat.workflow.firetasks.VaspTask` uses the fireworker `env` to know the VASP command via `vasp_cmd`. 

3. `fworker/example_job_template.yaml` is an example job template, which is adjusted by the queue adapter before being submitted to the queue when using `pybat qlaunch`. First, the job script has to load the python environment, so it can use `rlaunch` and knows the details of the Fireworks/FireTasks. Next, you need to load the modules that give you access to the VASP installation, if necessary. The jobscript then uses `sed` to adjust the number of the fireworker, nodes based on the number of nodes of the job. Finally, the script executes the `rlaunch` command, which communicates with the mongoDB database detailed in the launchpad and runs the Fireworks, if any.<br><br> `leibniz_job_template.sh` is the job template for the `leibniz` cluster of CalcUA, the organisation that manages the clusters of the university of Antwerp. It's mainly here so our students can copy it and adjust if necessary. 

4. `fworker/example_qadapter.yaml` is the queue adapter file, which details the variables which are plugged into the `$${variable}` locations in the job template. It only requires that you add a logdir to make sure all the fireworks output files are put in it. The other variables that the queue_adapter changes in the job_template are adjusted when using `pybat qlaunch`.

In order to get these files and start adjusting them, just open a terminal and execute the following:

```
wget https://www.dropbox.com/s/ruv8j51z8euy2nb/config.tar.gz
tar -xvzf config.tar.gz
rm config.tar.gz
```

Once you've made the necessary changes, use `pybat config` to configure pybat based on the configuration files. Be sure to check the help of each configure command with `-h` before using them.


