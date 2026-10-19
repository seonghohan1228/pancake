"""Check exact selection, paired branches, failure gating and exclusive outputs."""
import json
import os
from pathlib import Path
import runpy
import sys
from tempfile import TemporaryDirectory
import time
from unittest.mock import patch

root = Path(__file__).resolve().parents[1]
launcher = runpy.run_path(str(root/'cases/02f_lge_mixture_gas/launch_campaign.py'))
sleep = time.sleep

with TemporaryDirectory(dir=root/'tmp') as temporary:
    for fail in ('', 'P01'):
        campaign = Path(temporary)/('failed' if fail else 'complete')
        worker = campaign/'source/cases/02f_lge_mixture_gas/run_gas.py'
        worker.parent.mkdir(parents=True)
        worker.write_text("import argparse,os,sys\np=argparse.ArgumentParser();p.add_argument('--case');a,_=p.parse_known_args()\nsys.exit(3 if os.environ.get('PANCAKE_TEST_FAIL_CASE')==a.case else 0)\n")
        (campaign/'inputs').mkdir()
        selection = [f'P{i:02d}' for i in range(1,11)]
        experiment = dict(cases=[dict(case_id=key) for key in selection], pilot_case_ids=['P01','P07','P10'],
                          numerical_check_cases=dict(refinement=['P01','P10'],refinement_layers32=['P01']))
        (campaign/'inputs/experiment.json').write_text(json.dumps(experiment))
        with patch.object(sys,'argv',['launch','--root',str(campaign)]), \
             patch.object(time,'sleep',lambda _:sleep(.01)), \
             patch.dict(os.environ,{'PANCAKE_TEST_FAIL_CASE':fail}):
            assert launcher['main']() == (1 if fail else 0)
            saved=(campaign/'execution.json').read_bytes()
            state=json.loads(saved)
            jobs=list(state['jobs'].values())
            assert len(jobs) == (6 if fail else 26)
            assert all(j['case_id'] in selection for j in jobs)
            assert {j['case_id'] for j in jobs} == (set(experiment['pilot_case_ids']) if fail else set(selection))
            assert state['status'] == ('failed_needs_diagnosis' if fail else 'complete')
            assert all(j['state']!='running' for j in jobs)
            assert set(j['slot'] for j in jobs) <= set(range(6))
            for i,job in enumerate(jobs):
                assert job['branch'] in ('gas_on','gas_off')
                for other in jobs[i+1:]:
                    if job['slot']==other['slot']:
                        assert job['start_time']+job['elapsed_seconds'] <= other['start_time'] or other['start_time']+other['elapsed_seconds'] <= job['start_time']
            try:launcher['main']()
            except FileExistsError:pass
            else:raise AssertionError('Existing campaign must not be overwritten')
            assert (campaign/'execution.json').read_bytes()==saved
print('PASS: exact ten-case paired launcher, six checks, failure gating and preservation')
