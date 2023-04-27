#!/usr/bin/env python
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--nojec", action="store_true", help="whether or not to use the fit config with the JEC applied to msd."
)
parser.add_argument("--sample", choices=["WJets", "TTBar", "Combined"], default="WJets")
parser.add_argument("--years", nargs="+", default=["UL16preVFP", "UL16postVFP", "UL17", "UL18"])
parser.add_argument("--clear", action="store_true")
parser.add_argument("--prefix", "-p", default="python jetmass.py")
parser.add_argument("--pattern", "-P", default="")
parser.add_argument("--window", "-w", default=-1, type=int)
parser.add_argument("--exe", action="store_true")

args = parser.parse_args()

if "TMUX" not in os.environ:
    raise RuntimeError("This should only be run inside a tmux session!!")

window_nyears_panes = [
    k
    for k, v in [
            [int(f.replace("'", "")) for f in line.split(" ")]
            for line in subprocess.check_output(
                    ["tmux", "list-windows", "-F", "'#{window_index} #{window_panes}'"]
            ).splitlines()
        ]
    if v >= len(args.years)
]

if len(window_nyears_panes) == 0:
    raise RuntimeError("There is no tmux window with enough panes to fit all years!")
if args.window == -1:
    window_index = window_nyears_panes[0]
else:
    window_index = args.window
    if window_index not in window_nyears_panes:
        raise RuntimeError(
            "The provided window {} does not exist or has not enough panes to fit all years!".format(window_index)
        )

# cmd = "\"{} configs/".format(args.prefix) + args.sample + "_{YEAR}" + ("_noJEC" if args.nojec else "") + ".py\""
cmd = (
    '"{} '.format(args.prefix)
    + (
        args.pattern
        if args.pattern != ""
        else ("configs/" + args.sample + "_{YEAR}" + ("_noJEC" if args.nojec else "") + ".py")
    )
    + '"'
)

for iyear, year in enumerate(args.years):
    if args.clear:
        os.system("tmux send-keys -t :{}.{} \"C-c\"".format(window_index, iyear))
        os.system("tmux send-keys -t :{}.{} \"clear\" C-m".format(window_index, iyear))
        continue

    os.system("tmux send-keys -t :{}.{} {}".format(window_index, iyear, cmd.format(YEAR=year)))

    if args.exe:
        os.system("tmux send-keys -t :{}.{} \"C-m\"".format(window_index, iyear))
        # continue
