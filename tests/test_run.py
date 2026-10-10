"""Runner integration, resource telemetry and evidence preservation."""
import json
import numpy as np
import pytest
from prometheus.run import peak_memory_bytes, run_benchmark


def test_memory_telemetry():
    value=peak_memory_bytes()
    assert value is None or (isinstance(value,int) and value>0)


def test_pilot_archives_real_solves_and_refuses_overwrite(tmp_path):
    output=tmp_path/'pilot'
    summary=run_benchmark(output,repeats=1,meshes=(2.,),max_evaluations=1)
    assert summary['optimization'][0]['status']=='budget_exhausted'
    assert summary['optimization'][0]['evaluations']==1
    assert not any(v['passed'] for v in summary['verification'].values())
    assert not summary['failures']
    saved=json.loads((output/'summary.json').read_text())
    assert saved==summary
    assert (output/'manifest.json').is_file()
    assert np.load(output/'optimization_1'/'physical.npy').shape==(1600,)
    cost=summary['cost']
    assert cost['completed_factorizations']==1+cost['evaluation_factorizations']+2
    assert cost['completed_rhs_solves']==1+cost['evaluation_rhs_solves']+4
    before=(output/'summary.json').read_bytes()
    with pytest.raises(FileExistsError):
        run_benchmark(output,repeats=1,max_evaluations=1)
    assert (output/'summary.json').read_bytes()==before
