"""Strict resume validation with an explicitly reviewed code-transition bridge.

Resources are launch parameters, not analysis identity. All other identity fields
must match. Only exact old/new implementation fingerprints in the bridge permit
a code transition; unknown changes never inherit the old task-cache fingerprint.
"""


def resume_identity(previous, current, bridge):
    changed = {key for key in previous.keys() | current.keys()
               if previous.get(key) != current.get(key)}
    if not changed:
        return previous
    if changed == {'code_fingerprint'}:
        if (current['code_fingerprint'] == bridge.get('target_code_fingerprint')
                and previous['code_fingerprint'] in bridge.get('source_code_fingerprints', [])):
            return previous
    raise ValueError('Run inputs/code changed. Refusing incompatible resume; changed fields: '
                     + ', '.join(sorted(changed)))
