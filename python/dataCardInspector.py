#!/usr/bin/env coffea_python.sh
from simple_term_menu import TerminalMenu


class DataCardParser(object):
    def __init__(self, datacard_path):
        self.datacard_path = datacard_path
        self.clean_datacard = self.read_datacard()
        self.datacard = self.parse_datacard()

    def read_datacard(self):
        _clean_datacard = []
        with open(self.datacard_path, "r") as datacard_file:
            for line in datacard_file:
                clean = False
                l_ = list(line)
                while not clean:
                    for ichar in range(len(l_)):
                        if ichar == len(l_)-1:
                            clean = True
                        if l_[ichar] == " ":
                            if ichar < len(l_) and l_[ichar + 1] == " ":
                                l_.pop(ichar)
                                break
                        else:
                            continue
                _clean_datacard.append("".join(l_).strip())
        return _clean_datacard

    def parse_datacard(self):
        _datacard_dict = {"process_lines": [], "extArg": [], "nuisances": {}}

        observation_done = False
        for line in self.clean_datacard:
            if line.strip().startswith("#"):
                continue
            if any(max_num in line for max_num in ["imax", "jmax", "kmax"]):
                _datacard_dict[[max_num for max_num in ["imax", "jmax", "kmax"] if max_num in line][0]] = int(
                    line.split(" ")[1]
                )
                continue
            if "shapes" in line:
                _datacard_dict["shapes"] = line
                continue
            if "observation" in line:
                _datacard_dict["observations"] = [float(ob) for ob in line.split(" ") if not ob.isalpha()]
                observation_done = True
                continue
            if "bin" in line:
                bins = [_bin for _bin in line.split() if _bin != "bin"]
                if observation_done:
                    _datacard_dict["process_bins"] = bins
                else:
                    _datacard_dict["observations_bins"] = bins
                continue
            if "process" in line:
                _datacard_dict["process_lines"].append([_proc for _proc in line.split(" ") if _proc != "process"])
                continue
            if "rate" in line:
                _datacard_dict["rate"] = [_rate for _rate in line.split(" ") if _rate != "rate"]
                continue
            if "extArg" in line:
                _datacard_dict["extArg"].append(line.split(" ")[0])
                continue
            nuisance_name, nuisance_type = line.split(" ")[0:2]
            _datacard_dict["nuisances"].update(
                {nuisance_name: {"type": nuisance_type, "process_effect": [_eff for _eff in line.split(" ")[2:]]}}
            )
        return _datacard_dict

    def print_process(self, iproc, _print=False):
        if isinstance(iproc, str):
            iproc = self.datacard["process_lines"][0].index(iproc)
        if iproc >= self.datacard["jmax"]:
            print("provided process index exceeds number of processes in datacard!")
            return
        print_string = ""
        print_string += "%s (%s)\n" % tuple(procline[iproc] for procline in self.datacard["process_lines"])
        print_string += "rate: %s\n" % self.datacard["rate"][iproc]
        print_string += "nuisances:\n"
        for nu_name in self.datacard["nuisances"]:
            nu_dict = self.datacard["nuisances"][nu_name]
            if not nu_dict["process_effect"][iproc].isalpha():
                print_string += "%s (%s) = %s\n" % (nu_name, nu_dict["type"], nu_dict["process_effect"][iproc])

        if _print:
            print(print_string)
        return print_string

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("datacard", help="Path to datacard.")

    args = parser.parse_args()

    datacard = DataCardParser(args.datacard)

    inspector_options = ["channels", "exit"]
    exit_ = False

    main_menu_cursor = "> "
    main_menu_cursor_style = ("fg_green", "bold")
    main_menu_style = ("bg_red", "fg_green")
    main_menu_exit = False
    
    inspector_menu = TerminalMenu(
        inspector_options,
        # multi_select=True,
        # show_multi_select_hint=True,
        cursor_index=inspector_options.index("exit"),
        menu_cursor=main_menu_cursor,
        menu_cursor_style=main_menu_cursor_style,
        menu_highlight_style=main_menu_style,
        cycle_cursor=True,
    )
    
    channel_options = datacard.datacard["process_lines"][0] + ["back"]
    channel_menu = TerminalMenu(
        channel_options,
        multi_select=True,
        show_multi_select_hint=True,
        # cursor_index=channel_options.index("back"),
        menu_cursor=main_menu_cursor,
        menu_cursor_style=main_menu_cursor_style,
        menu_highlight_style=main_menu_style,
        cycle_cursor=True,
        # preview_command=datacard.print_process,
        # preview_size=0.75
    )
    while not exit_:
        menu_choice = inspector_options[inspector_menu.show()]
        # for job in menu_choices:
        if menu_choice == "channels":
            channel_menu_choices_idxs = channel_menu.show()
            channel_menu_choices = list(channel_menu.chosen_menu_entries)
            for channel_job in channel_menu_choices:
                if channel_job in "back":
                    break
                datacard.print_process(channel_job, _print=True)

        if menu_choice == "exit":
            exit_ = True
                
