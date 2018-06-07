# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob

from fireworks import Firework, LaunchPad, PyTask, FWorker, \
    Workflow
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

"""
Workflow setup for the pybat package.

"""

# Set up the Launchpad for the workflows
LAUNCHPAD = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                      username="mbercx", password="quotastests")