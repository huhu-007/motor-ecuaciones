import threading
import time
import unittest

from runner import Busy, Runner


class TestRunner(unittest.TestCase):
    def setUp(self):
        self.r = Runner(max_concurrent=1, queue_wait=0.3)

    def test_solve_ok(self):
        out = self.r.run("solve", "x^2 - 2 = 0", 20, timeout=20)
        self.assertEqual(out["status"], "solved")
        self.assertEqual(out["roots"][0]["exactness"], "exact_proven")

    def test_invalid_input(self):
        out = self.r.run("solve", "x^^^2 ===", 20, timeout=20)
        self.assertEqual(out["status"], "invalid_input")

    def test_analyze(self):
        out = self.r.run("analyze", "√3+√2", 60, timeout=20)
        self.assertEqual(out["exactness"], "exact_proven")

    def test_timeout_kills_process(self):
        t = time.time()
        out = self.r.run("solve", "sin(x)=0", 30, timeout=0.05)
        self.assertEqual(out["status"], "timeout")
        self.assertLess(time.time() - t, 3)

    def test_busy(self):
        res = {}
        th = threading.Thread(target=lambda: res.update(self.r.run("solve", "sin(x)=0", 30, timeout=20)))
        th.start()
        time.sleep(0.1)
        with self.assertRaises(Busy):
            self.r.run("solve", "x-1=0", 10, timeout=5)
        th.join()
        self.assertEqual(res["status"], "solved")


if __name__ == "__main__":
    unittest.main()
