from multiprocessing import Process, Queue


class Multiprocessor:

    def __init__(self):
        self.processes = []
        self.procid = []
        self.queue = Queue()

    @staticmethod
    def _wrapper(func,procID, queue, args, kwargs):
        ret = func(*args, **kwargs)
        queue.put([ret,procID])

    def run(self, func, procID, *args, **kwargs):
        args2 = [func, procID, self.queue, args, kwargs]
        p = Process(target=self._wrapper, args=args2)
        self.processes.append(p)
        self.procid.append(procID)
        p.start()

    def wait(self):
        rets = [None] * len(self.processes)
        for p in self.processes:
            _ret = self.queue.get()
            rets[_ret[1]] = _ret[0]
        for p in self.processes:
            p.join()
        return rets
