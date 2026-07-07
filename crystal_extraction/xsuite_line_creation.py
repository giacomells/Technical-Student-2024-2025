"""Tool to retrieve lattices from gitlab and create xsuite lines."""

import typing as t

import requests
import xtrack as xt
from cpymad.madx import Madx

DEFAULT_URLS = [
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/SPS_LS2_2020-05-26.seq",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/strengths/lhc_q20.str",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/toolkit/macro.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_classes.madx",
    "https://gitlab.cern.ch/acc-models/acc-models-sps/-/raw/2021/aperture/aperturedb_elements.madx",
]

def create_xsuite_line(sequence: str = 'sps', file_path: t.Optional[str] = None, urls: t.Optional[list] = None,
                        mad_inputs: t.Optional[list] = None, mad_stdout = None, from_madx_sequence_kwargs: t.Optional[dict] = None) -> xt.Line:
    """
    Create an xsuite line from a MAD-X sequence.

    Parameters:
    sequence (str): The name of the MAD-X sequence to use. Default is 'sps'.
    file_path (str): The file path to save the xsuite line as a JSON file. Default is None.
    urls (list): A list of URLs to retrieve MAD-X input files. Default is DEFAULT_URLS.
    mad_inputs (list): Additional MAD-X input commands. Default is None.
    mad_stdout: The standard output for MAD-X. Default is None.

    Returns:
    xt.Line: The created xsuite line.
    """
    if urls is None:
        urls = DEFAULT_URLS
    if from_madx_sequence_kwargs is None:
        from_madx_sequence_kwargs = {"install_apertures": True, "deferred_expressions": True}

    mad_ = Madx(stdout=mad_stdout)
    mad_.beam()

    for url_ in urls:
        response = requests.get(url_)
        response.raise_for_status()
        mad_.input(response.text)

    if mad_inputs is not None:
        for input_ in mad_inputs:
            mad_.input(input_)

    mad_.use(sequence)
    xline = xt.Line.from_madx_sequence(mad_.sequence[sequence],
                                       **from_madx_sequence_kwargs)

    if file_path is not None:
        xline.to_json(file_path)

    return xline

