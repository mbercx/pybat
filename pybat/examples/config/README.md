# Configuring pybat for workflows via fireworks

In this directory you can find some example configuration files for the pybat workflow setup. These files correspond to the typical [Fireworks](https://materialsproject.github.io/fireworks/index.html) YAML files that describe both the fireworker and the launchpad. Here's an overview:

1. `fworker/example_fworker.yaml`: This file configures the fireworker, i.e. the cluster that performs the calculations. Important here is that pybat uses the fireworker category to make sure that Fireworks that require a specific amount of nodes are only picked up by the corres job. Also: the `pybat.workflow.firetasks.VaspTask` uses the fireworker `env` to know the VASP command via `vasp_cmd`.

2. `fworker/example_job_template.yaml`: This is an example job template for the `leibniz` cluster of CalcUA, the organisation that manages the clusters of the university of Antwerp. It's mainly here so our students can copy it and adjust. Important here is that the jobscript uses sed to adjust the number of the fireworker, nodes based on the number of nodes of the job.

3. `fworker/example_job_template.yaml`: This file only requires that you add a logdir to make sure all the fireworks output files are put in it. The other variables that the queue_adapter changes in the job_template are adjusted when using `pybat qlaunch`.

4. `launchpad/example_launchpad.yaml`: This is an example of the launchpad file, i.e. the file that contains all the information about your mongoDB database which is used to store the workflows. 

In order to get these files and start adjusting them, just open a terminal and execute the following:

```
wget https://www.dropbox.com/s/ruv8j51z8euy2nb/config.tar.gz
tar -xvzf config.tar.gz
rm config.tar.gz
```

Once you've made the necessary changes