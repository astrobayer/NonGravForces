"""Is every ADMISSIBLE OPTION VALUE described where a user would look for it?

Options.validateValue holds the authority for what each dotted field accepts. This asks two
things of every value in it:
  1. is it named in Options.m's own help, so `help NonGravForces.Options` lists it; and
  2. is it named in the help of the UNIT that consumes it, which is where its meaning lives
     (drag.accommodation "sesam" is explained by kernel.accommodation, not by Options).
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, '+NonGravForces')

# which unit's help is expected to explain each family's values
CONSUMER = {
    'drag_density': '+env/density.m',
    'drag_wind': '+env/wind.m',
    'drag_gsi': '+kernel/gsiCoefficients.m',
    'drag_accommodation': '+kernel/accommodation.m',
    'srp_eclipse': '+env/eclipseFactor.m',
    'erp_backend': '+env/earthRadiance.m',
    'thermal_model': '+geom/thermalState.m',
    'geom_selfShadow': '+geom/selfShadow.m',
    'sun_source': '+env/radiationSources.m',
    'partials_method': 'Model.m',            # Model.partials selects it; kernel.partials is the analytic form itself
}


def main():
    src = open(os.path.join(PKG, 'Options.m'), encoding='utf-8').read()
    block = re.search(r"sets\s*=\s*struct\((.*?)\);", src, re.S)
    if not block:
        print('could not find Options.validateValue sets')
        return 2
    pairs = re.findall(r"'(\w+)',\s*\[([^\]]*)\]", block.group(1))
    helptext = '\n'.join(l for l in src.splitlines() if l.lstrip().startswith('%'))

    findings = 0
    for key, raw in pairs:
        values = re.findall(r'"([^"]+)"', raw)
        missing_here = [v for v in values if not re.search(r'"%s"|\b%s\b' % (re.escape(v), re.escape(v)), helptext)]
        if missing_here:
            findings += 1
            print(f'Options.m help does not name {key} value(s): {", ".join(missing_here)}')
        rel = CONSUMER.get(key)
        if not rel:
            continue
        path = os.path.join(PKG, *rel.split('/'))
        if not os.path.isfile(path):
            print(f'  (consumer {rel} for {key} not found)')
            continue
        csrc = open(path, encoding='utf-8').read()
        chelp = '\n'.join(l for l in csrc.splitlines() if l.lstrip().startswith('%'))
        # "none" means the consuming unit is NOT CALLED at all, so its own help has no reason to
        # name it -- erp.backend "none" is handled by radiationSources, not by earthRadiance.
        missing_there = [v for v in values
                         if v != 'none' and not re.search(r'\b%s\b' % re.escape(v), chelp)]
        if missing_there:
            findings += 1
            print(f'{rel} help does not explain {key} value(s): {", ".join(missing_there)}')

    print(f'\n{findings} finding(s) over {len(pairs)} option families')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
