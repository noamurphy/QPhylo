"""Paste into Qollab's main file. Run one selected check at a time."""

import json
from math import pi
import time
import qiskit
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.providers.jobstatus import JobStatus

test_name = "bit_order_q0"
# Next: bit_order_q2, rx_pi, rz_phase, zz_phase.
shots = 1000
backend_label = "UNRECORDED: copy the exact simulator/noise selection here"

theta = Parameter("theta")
if test_name in ("bit_order_q0", "bit_order_q2"):
    qc = QuantumCircuit(3, 3)
    qubit = 0 if test_name == "bit_order_q0" else 2
    qc.x(qubit)
    expected = "001" if qubit == 0 else "100"
elif test_name == "rx_pi":
    qc = QuantumCircuit(1, 1)
    qc.rx(theta, 0)
    qc = qc.assign_parameters({theta: pi})
    expected = "1"
elif test_name == "rz_phase":
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    qc.rz(theta, 0)
    qc.h(0)
    qc = qc.assign_parameters({theta: pi})
    expected = "1"
elif test_name == "zz_phase":
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.rz(theta, 1)
    qc.cx(0, 1)
    qc.h(0)
    qc = qc.assign_parameters({theta: pi})
    expected = "01"
else:
    raise ValueError(f"Unknown test: {test_name}")
qc.measure(list(range(qc.num_qubits)), list(range(qc.num_qubits)))

print("Check:", test_name)
print("Qiskit version:", qiskit.__version__)
print("Backend type:", type(backend).__name__)
print("Backend selection:", backend_label)
print("Expected ideal bitstring:", expected)
print(qc)
started = time.monotonic()
job = backend.run(qc, shots=shots)
print("Job ID:", job.job_id())
last_status = None
while True:
    status = job.status()
    if status != last_status:
        print("Job status:", status)
        last_status = status
    if status in (JobStatus.DONE, JobStatus.ERROR, JobStatus.CANCELLED):
        break
    if time.monotonic() - started >= 600:
        raise TimeoutError("Stopped waiting; the remote job may still be running.")
    time.sleep(5)
if status != JobStatus.DONE:
    raise RuntimeError(f"Job failed to complete: {status}")
counts = job.result().get_counts()
if sum(counts.values()) != shots:
    raise ValueError("Returned measurement count does not match requested shots.")

report = {"test": test_name, "shots": shots, "expected": expected,
          "counts": counts, "expected_fraction": counts.get(expected, 0) / shots,
          "backend_label": backend_label, "backend_type": type(backend).__name__,
          "qiskit_version": qiskit.__version__, "job_id": str(job.job_id()),
          "elapsed_seconds": time.monotonic() - started}
print("Counts:", counts)
print(f"Expected outcome observed: {report['expected_fraction']:.1%}")
print("REPORT_JSON:", json.dumps(report, sort_keys=True))
print("PASS: job completed and returned all shots. Compare outcome with expectation.")
