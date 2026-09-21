import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'launch_support'))
sys.path.insert(0, str(ROOT/'scripts'))
from resource_resume import resume_identity
from common import digest, sha, tree_hash


class ResourceResume(unittest.TestCase):
    def setUp(self):
        self.bridge = json.loads((ROOT/'launch_support/resource_resume_bridge.json').read_text())
        self.old = dict(code_fingerprint=self.bridge['source_code_fingerprints'][0],
                        selection_sha256='selection', settings_sha256='settings',
                        conda_packages_sha256='environment', work_dir='same-work')
        self.new = dict(self.old, code_fingerprint=self.bridge['target_code_fingerprint'])

    def test_reviewed_transition_preserves_cache_identity(self):
        self.assertEqual(resume_identity(self.old, self.new, self.bridge), self.old)

    def test_identical_resume(self):
        self.assertEqual(resume_identity(self.new, self.new, self.bridge), self.new)

    def test_changed_inputs_and_environment_rejected(self):
        for field in ('selection_sha256', 'settings_sha256', 'conda_packages_sha256', 'work_dir'):
            with self.assertRaises(ValueError):
                resume_identity(self.old, dict(self.new, **{field: 'different'}), self.bridge)

    def test_unknown_code_rejected(self):
        with self.assertRaises(ValueError):
            resume_identity(self.old, dict(self.new, code_fingerprint='unreviewed'), self.bridge)

    def test_shared_tool_change_cannot_use_resource_only_bridge(self):
        actual = digest(dict(scripts=tree_hash(ROOT/'scripts'), workflow=sha(ROOT/'main.nf'),
                             config=sha(ROOT/'nextflow.config'), cluster=tree_hash(ROOT/'conf'),
                             launcher=sha(ROOT/'run_validation.py'), environment=sha(ROOT/'environment.yml'),
                             resume_helper=sha(ROOT/'launch_support/resource_resume.py')))
        self.assertNotEqual(actual, self.bridge['target_code_fingerprint'])
        with self.assertRaises(ValueError):
            resume_identity(self.old, dict(self.new, code_fingerprint=actual), self.bridge)


if __name__ == '__main__': unittest.main()
