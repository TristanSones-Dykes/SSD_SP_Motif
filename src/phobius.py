from bs4 import BeautifulSoup as bs
import mechanicalsoup
import pandas as pd
import os
import signalp as sp


# ------ Description ------
# This script contains functions to run phobius transmembrane
# prediction on a list of proteins in FASTA format by web scraping
# the form of the online tool

# url for phobius
url = "https://phobius.sbc.su.se/index.html"

# function to run phobius on a FASTA formatted list of proteins
# reads from file by default
def phobius(origin: str, isString: bool = False) -> pd.DataFrame:
    if isString:
        fileString = origin
    else:
        with open(origin, "r") as f:
            fileString = f.read()

    #instantiate browser, open url, select form
    browser = mechanicalsoup.StatefulBrowser(soup_config={'features': 'html.parser'})
    browser.open(url)
    browser.select_form('form')

    # fill in form and check
    browser["format"] = "nog"
    if isString:
        browser["protseq"] = fileString
    else:
        browser.form.set("protfile", origin)

    response = browser.submit_selected()

    # parse response to get all data in <pre> tag
    soup = bs(response.text, 'html.parser')
    out = soup.find('pre')

    # right -> left:
    # splits output text into list of protein outputs 
    # then splits each protein output into list of non-empty lines and collapses whitespace
    out_split = [[y.split() for y in x.split("\n") if y != ''] for x in out.text.split("//")][:-1]

    assert len(out_split) == len(fileString.split(">")) - 1, "number of proteins in output does not match number of proteins in input"

    # split into predicted and unpredicted
    predicted = []
    unpredicted = []
    for protein in out_split:
        if 'TRANSMEM' in [x[1] for x in protein]:
            predicted.append(protein)
        else:
            unpredicted.append(protein)

    # populate dataframe with predicted transmembrane window sizes
    window_df = pd.DataFrame(columns = ["ID", "start", "end"], dtype = int)
    for protein in predicted:
        # extract category column
        col_2 = [x[1] for x in protein]

        # extract protein ID, and start and end of transmembrane region
        proteinID = col_2[0]
        transmem_row = col_2.index("TRANSMEM")
        start = int(protein[transmem_row][2])
        end = int(protein[transmem_row][3])
        
        # append to dataframe
        new_row = pd.DataFrame([[proteinID, start, end]], columns = ["ID", "start", "end"])
        window_df = pd.concat([window_df, new_row], ignore_index=True)
    
    #delete temp file
    if isString:
        os.remove("temp.fasta")

    return window_df