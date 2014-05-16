# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from collections import OrderedDict
from uuid import uuid1
from time import time
from pprint import pformat
from datetime import datetime


class AutoVivification(dict):

    """Implementation of perl's autovivification feature. Checkout http://stackoverflow.com/a/652284/280182"""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class TimeMachine(object):

    """docstring for TimeMachine"""

    def __init__(self):
        super(TimeMachine, self).__init__()
        self.history = OrderedDict()

    def __call__(self, do=None, undo=None, bookmark=None):
        do()
        current_time = time()
        if bookmark is None:
            entry_id = uuid1()
        else:
            entry_id = bookmark
        # make sure that entry is added to the end of history
        self.history.pop(entry_id, None)
        self.history[entry_id] = {'unix_epoch':
                                  current_time, 'undo': undo, 'redo': do}
        return entry_id

    def __str__(self):
        info = '\n'
        for uuid, entry in self.history.iteritems():
            info += datetime.fromtimestamp(entry['unix_epoch']
                                           ).strftime('%Y-%m-%d %H:%M:%S') + '\n'
            info += 'undo: ' + \
                ' '.join([str(elem)
                         for elem in (entry['undo'].func, entry['undo'].args, entry['undo'].keywords)]) + '\n'
            info += 'redo: ' + \
                ' '.join([str(elem)
                         for elem in (entry['redo'].func, entry['undo'].args, entry['undo'].keywords)]) + '\n'
        return info

    def __repr__(self):
        return self.__str__()

    def undo(self, bookmark=None):
        if bookmark is None:
            (uuid, entry) = self.history.popitem()
            entry['undo']()
        elif bookmark in self.history.keys():
            uuid = False
            while uuid is not bookmark:
                (uuid, entry) = self.history.popitem()
                entry['undo']()
        else:
            raise Exception(
                'Provided bookmark %s cannot be found in the time machine.')

    def reset(self):
        self.undo(bookmark=self.history.keys()[0])


def partition(lst, n):
    """Partition a list into n bite size chunks."""
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]