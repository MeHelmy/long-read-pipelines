import os

# Data manipulation
import numpy as np
import pandas as pd
import random as rnd
import tabulate
import re
from enum import Enum
import glob


# Plotting
import matplotlib.pyplot as plt
import plotly.express as px
import folium
import plotly.graph_objects as go
import math


# HTML
from markupsafe import Markup
import jinja2
import io
from io import StringIO
import base64


# Argument Parsing
import argparse


''' 
matplotlib render settings
'''
# Make big figures:
gFIG_SIZE_in = [14, 10]

# Set plotting defaults:
gPLOT_PARAMS = {
    "legend.fontsize": "x-large",
    "figure.figsize": gFIG_SIZE_in,
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large"
}
plt.rcParams.update(gPLOT_PARAMS)

# Some single-place definitions of sizes for plots / figures:
gFONT_SIZE_UNITS = "pt"
gTITLE_FONT_SIZE = 36
gAXIS_LABEL_FONT_SIZE = 24
gTICK_LABEL_FONT_SIZE = 16
gTEXT_FONT_SIZE = 16

# To track global figure number creation:
gFIG_NUM = 0


def fix_plot_visuals(fig,
                     titlesize=gTITLE_FONT_SIZE,
                     labelsize=gAXIS_LABEL_FONT_SIZE,
                     ticklabelsize=gTICK_LABEL_FONT_SIZE,
                     textsize=gTEXT_FONT_SIZE,
                     tight_rect=None):
    """Fix the plot elements to be appropriate sizes for a slide / presentation."""

    if not textsize:
        textsize = ticklabelsize

    for ax in fig.get_axes():

        for ticklabel in (ax.get_xticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for ticklabel in (ax.get_yticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for c in ax.get_children():
            if c.__class__ == matplotlib.text.Text:
                c.set_fontsize(textsize)

        ax.xaxis.get_label().set_fontsize(labelsize)
        ax.yaxis.get_label().set_fontsize(labelsize)
        ax.title.set_fontsize(titlesize)

    for c in fig.get_children():
        if c.__class__ == matplotlib.legend.Legend:
            c.prop.set_size(ticklabelsize)
            c.get_title().set_size(ticklabelsize)

    if tight_rect:
        fig.tight_layout(rect=tight_rect)
    else:
        fig.tight_layout()
    
    if fig._suptitle:
        sup_title = fig._suptitle.get_text()
        fig.suptitle(sup_title, fontsize=titlesize)
    
    # Make it so we can actually see what's happening on the plots with "dark mode":
    fig.patch.set_facecolor("white")


'''
Coverage Plot
'''
# Function to go through each .bed.gz file, iterate through it and bin it, and plot it
# to a larger table which is plotted later.
def plot_coverage(directory, sample_name, bin_width=500):
    
    # only looking for .bed.gz files
    ext = (".bed.gz")

    # Set up DF to hold all data
    beds = pd.DataFrame(columns=["stop", "depth"])

    # Reading and sorting files
    file_list = os.listdir(directory)
    dir_len = len(file_list)
    print(f"Files to be plotted: {file_list}")
    
    # Plot setup
    fig, ax = plt.subplots(1,1, figsize=(15, 9))
    plt.title(f"Coverage Plot of Sample {sample_name}", pad=12, fontsize=12)
    plt.xlabel("Contig (bp)", labelpad=10, fontsize=12)
    plt.ylabel("Depth", labelpad=10, fontsize=12)
    color = plt.cm.viridis(np.linspace(0, 1, dir_len))
    tick_positions = []
    contigs = []
    bin_max = 0
 
    # Iterate over each bed.gz in test_data folder
    for idx, file in enumerate(file_list):

        if file.endswith(ext):
            # Reading current .bed.gz file
            f = os.path.join(directory, file)
            #print(f"\n{idx} Current file: {f}")

            # Create DataFrame per .bed.gz
            bed_df = pd.read_csv(f, delimiter="\\t", names=["contig", "start", "stop", "depth"], engine="python")
            
            # Get bins
            #print("Getting bins...")
            start_min = bed_df["start"].min()
            stop_max = bed_df["stop"].max()
            bins = np.arange(start_min, stop_max+1, bin_width)
            values = np.zeros(len(bins)) 
            #print(f"Max Pos: {stop_max}")

            # Iterrate through each DataFrame and populate bins
            #print("Populating bins...")
            for _, row in bed_df.iterrows():
                avg = (row["stop"]+row["start"])/2
                index = int(np.floor(avg/bin_width))
                values[index]+=row["depth"]

            # Append new data to DF
            #print("Plotting data...")
            if(idx == 0):
                ax.plot(bins, values, ".", c=color[idx])
                tick_positions.append(math.floor(stop_max/2)+bin_max)
                bin_max = max(bed_df["stop"])
            else:
                ax.plot(bins+bin_max, values, ".", c=color[idx])
                tick_positions.append(math.floor(stop_max/2)+bin_max)
                bin_max = bin_max + max(bed_df["stop"])                      
                               
            # Saving xtick data
            
            contigs.append(bed_df["contig"][0])
            
            print(f"{idx+1}/{dir_len} files completed!")
                
    # Setting xticks
    ax.set_xticks(ticks=tick_positions)
    ax.set_xticklabels(labels=contigs, rotation=90)
    
    fix_plot_visuals(fig)
            
    return fig

# Function to convert given plot to b64 using I/O buffer
def img_to_b64(plot):
    
    # Create I/O buffer to save image in bytes
    stringIObytes = io.BytesIO()
    plot.savefig(stringIObytes, format="jpeg")
    
    # Retrieve byte-string and encode it in base64
    stringIObytes.seek(0)
    plot_b64 = base64.b64encode(stringIObytes.read())
    
    return plot_b64


'''
Drug Resistance Table
'''
def create_drug_table(file):
    if not file:
        resistances = ["UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED"]
    else:
        data = open(file, 'r').read()
        resistances = list(get_drug_resistance(data, None, None, do_print=True))
    
    resistances_tbl = pd.DataFrame(columns = ["Chloroquine", "Pyrimethamine", "Sulfadoxine", "Mefloquine", "Artemisinin", "Piperaquine"])
    resistances = map(str, resistances)
    resistances = [s.replace('DrugSensitivity.','') for s in list(resistances)]
    
    return resistances

# Set up some functions to determine drug sensitivity:
# vars for markers in drug_resistance_report.txt
present = '+'
absent = '-'

pchange_regex = re.compile(r"""p\.([A-z]+)([0-9]+)([A-z]+)""")
AA_3_2 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

DrugSensitivity = Enum("DrugSensitivity", ["UNDETERMINED", "SENSITIVE", "RESISTANT"])

def parse_pchange(pchange):
    old, pos, new = pchange_regex.match(pchange).groups()
    return AA_3_2[old.upper()], int(pos), AA_3_2[new.upper()]

def get_chloroquine_sensitivity(dr_report):
# Locus utilized: PF3D7_0709000 (crt)
# Codon: 76
# Workflow:
# Step Genetic change Interpretation Classification
# 1 76 K/T heterozygote Heterozygous mutant Undetermined
# 2 76 missing Missing Undetermined
# 3 K76 Wild type Sensitive
# 4 76T Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfcrt PF3D7_0709000"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 76:
                if (old == "L") and (new == "T") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (new == "T") and (marker == present):
                    return DrugSensitivity.RESISTANT
                
    return DrugSensitivity.UNDETERMINED
                
def get_pyrimethamine_sensitivity(dr_report):
# Locus utilized: PF3D7_0417200 (dhfr)
# Codon: 108
# Workflow:
# Step Genetic change Interpretation Classification
# 1 108 S/N heterozygote Heterozygous mutant Undetermined
# 2 108 missing Missing Undetermined
# 3 S108 Wild type Sensitive
# 4 108N Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfdhfr PF3D7_0417200"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 108:
                if (old == "S") and (new == "N") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (old == "S") and (new == "N") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "N") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif marker == absent:
                    return DrugSensitivity.SENSITIVE
                
    return DrugSensitivity.UNDETERMINED


def get_sulfadoxine_sensitivity(dr_report):
# Locus utilized: PF3D7_0810800 (dhps)
# Codon: 437
# Workflow:
# Step Genetic change Interpretation Classification
# 1 437 A/G heterozygote Heterozygous mutant Undetermined
# 2 437 missing Missing Undetermined
# 3 A437 Wild type Sensitive
# 4 437G Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())

        # We only care about this locus for this drug:
        if line.startswith("pfdhps PF3D7_0810800"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 437:
                if (old == "A") and (new == "G") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (old == "A") and (new == "G") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "G") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif marker == absent:
                    return DrugSensitivity.SENSITIVE
                
    return DrugSensitivity.UNDETERMINED

def get_mefloquine_sensitivity(dr_report):
# Locus utilized: PF3D7_0523000 (mdr1)
# Codons: Amplification status of whole gene
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Missing Missing Undetermined
# 2 Heterozygous duplication Heterozygous mutant Undetermined
# 3 Single copy Wild type Sensitive
# 4 Multiple copies Mutant Resistant

    # Currently we can't determine this.
    # We need to get CNV calling working first.
                
    return DrugSensitivity.UNDETERMINED

def get_artemisinin_sensitivity(dr_report):
# Locus utilized: PF3D7_1343700 (kelch13)
# Codons: 349-726 (BTB/POZ and propeller domains)
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Homozygous non-synonymous mutations in the kelch13 BTB/POZ and propeller
# domain classified by the World Health Organisation as associated with delayed
# parasite clearance 
# Mutant – associated with delayed clearance Resistant
# 2 Heterozygous non-synonymous mutations in the kelch13 BTB/POZ and
# propeller domain classified by the World Health Organisation as associated
# with delayed parasite clearance
# Mutant - heterozygous Undetermined
# 3 578S as homozygous Mutant - not associated Sensitive
# 4 Any missing call in amino acids 349-726 Missing Undetermined
# 5 No non-reference calls in amino acids 349-726 Wild-type Sensitive
# 6 otherwise Mutant - not in WHO list Undetermined

    has_variants = False

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfkelch13 PF3D7_1343700"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            
            has_non_ref = False
            has_variants = False
            if 349 <= pos <= 726:
                if (old != new) and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "S") and (marker == present):
                    return DrugSensitivity.SENSITIVE
                elif marker == present:
                    has_non_ref = True
                has_variants = True
    if (has_variants) and (not has_non_ref):
        return DrugSensitivity.SENSITIVE
    
    return DrugSensitivity.UNDETERMINED

def get_piperaquine_sensitivity(dr_report):
# Loci utilized: PF3D7_1408000 (plasmepsin 2) and PF3D7_1408100 (plasmepsin 3)
# Codons: Amplification status of both genes
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Missing Missing Undetermined
# 2 Heterozygous duplication Heterozygous mutant Undetermined
# 3 Single copy Wild type Sensitive
# 4 Multiple copies Mutant Resistant

    # Currently we can't determine this.
    # We need to get CNV calling working first.
                
    return DrugSensitivity.UNDETERMINED

def get_drug_resistance(data, sample_id, sample_df, do_print=False):
    # Get the GSURI to the drug resistance report:
    #row = reportable_samples.loc[reportable_samples["entity:sample_id"] == sample_id]
    #dr_report_gs_url = row.iloc[0]["drug_res_report"]

    #blob = storage.Blob.from_string(dr_report_gs_url)
    
    # Download the file contents:
    #dr_report_contents = blob.download_as_text(client=my_storage_client)
    
    chloroquine = get_chloroquine_sensitivity(data)
    pyrimethamine = get_pyrimethamine_sensitivity(data)
    sulfadoxine = get_sulfadoxine_sensitivity(data)
    mefloquine = get_mefloquine_sensitivity(data)
    artemisinin = get_artemisinin_sensitivity(data)
    piperaquine = get_piperaquine_sensitivity(data)

   # if do_print:
   #     print("\t".join([
   #         "Sample\t", "Chloroquine", "Pyrimethamine", "Sulfadoxine", "Mefloquine", "Artemisinin", "Piperaquine"
   #     ]))
   #     print(f"{sample_id}\t", end='')
   #     for r in [chloroquine, pyrimethamine, sulfadoxine, mefloquine, artemisinin, piperaquine]:
   #         print(f"{r.name}\t", end='')
   #     print()
        
    return chloroquine, pyrimethamine, sulfadoxine, mefloquine, artemisinin, piperaquine
 
 
'''
Map
'''
def create_map(coordinates, sample_name):
    m = folium.Map(location=coordinates, zoom_start = 5)

    # Check if make_default == True -- if so, do not make a marker
    #if not make_default:
    #    folium.Marker(location=coordinates, popup = ('Sample: '+sample_name), icon=folium.Icon(color='red',prefix='fa',icon='circle'), parse_html=True).add_to(m)

    if coordinates != [0,0]:
        folium.Marker(location=coordinates, popup = ('Sample: '+sample_name), icon=folium.Icon(color='red',prefix='fa',icon='circle'), parse_html=True).add_to(m)

    m.get_root().width = "473px"
    m.get_root().height = "397px"
    map_html = m.get_root()._repr_html_()
    
    #with open('/templates/map.html', mode='w+') as f:
    #    f.write(map_html)
    
    return map_html


'''
Q-Score Plot
'''
def create_qscore_plot():
    '''
    Function to create a q-score index vs. reads plot
    
    Plot is returned as an HTML file. Currently, this plot is written
    in Javascript. This function is unused in report generation.
    '''
    
    x=[5, 7, 10, 12, 15] # index
    y=[3413205, 3413204, 3413073, 3062945, 2120402] # q scores
    
    df = pd.DataFrame({
        'Q Score': x,
        'Number of Reads': y
    })
    
    fig1 = px.line(df,x="Q Score",y="Number of Reads", width=450, height=400)
    
    outpath="github/lr-malaria-automation-report/report_files/templates/qscore_plot.html"
    return fig1.write_html(outpath,
                full_html=False,
                include_plotlyjs=True)


'''
Report Generation & Templating
'''
class Sample:
    '''
    This class holds all variables that will be used on the summary page of the report.
    
    Additionally, this class holds variables that are used on both the summary and analysis pages,
    such as the sample name.
    '''
    
    def __init__(self, sample_name, hrp2, hrp3, qc_status, drug_res, info, _map, location_info):
        '''This function defines the class variables and retrieves them from their respective functions.'''
        
        self.sample_name = sample_name
        self.hrp2 = hrp2
        self.hrp3 = hrp3
        self.qc_status = qc_status
        self.drug_res = drug_res
        self.info = info
        self.map = _map
        self.location_info = location_info

class Analysis:
    '''
    This class holds all variables used on the analysis page of the report.
    
    Additionally, it passes in variables needed for plot generation in Javascript.
    '''
    
    def __init__(self, sequencing_summary, qscore_reads, qscore_scores, active_channels, coverage_plot):
        '''This function defines the variables used and retrieves them from their respective functions.'''
        
        self.sequencing_summary = sequencing_summary
        self.scores = qscore_scores
        self.reads = qscore_reads
        self.active_channels = active_channels
        self.coverage_plot = coverage_plot
        # self.reference_info = reference_info


def create_report(sample, analysis):
    '''
    This function handles the final steps of report generation.
    
    It pools all variables and data needed for the report and creates two objects:
    analysis and summary. These objects hold the variables and are passed to their respective
    page templates (analysis or summary).
    '''
    print("path: "+os.path.dirname(os.path.realpath(__file__)))
    print("cwd: "+os.getcwd())

    # creating summary page
    templateLoader = jinja2.FileSystemLoader(searchpath='/report-files/templates/')
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = 'report.html' # may need to change if file is moved
    template = templateEnv.get_template(TEMPLATE_FILE)
    output = template.render(sample=sample, analysis=analysis)

    print(os.listdir())
    file_name = os.path.join(os.getcwd(),sample.sample_name+'_lrma_report.html')
    print(file_name)

    with open(file_name, 'w') as fh:
        fh.write(output)  

    print(os.listdir())
        
    print('Report generated!')

if __name__ == '__main__':
    '''Reports will be generated in the current working directory.
    
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Set up argument parsing
    '''
    # define parser object
    parser = argparse.ArgumentParser()

    # define accepted input arguments
    ''' Summary Page '''
    # required inputs
    
    # Sample Info
    parser.add_argument("--sample_name", help="name of sequenced sample", required=True)
    parser.add_argument("--upload_date", help="date sample was sequenced and uploaded", nargs='+', required=True)
    parser.add_argument("--species", help="species of sample", nargs='+', default="P. falciparum")
    parser.add_argument("--aligned_coverage", help="number of times the bases in the sequenced reads cover the target genome", required=True, type=float) # check -- fold coverage
    parser.add_argument("--aligned_read_length_n50", help="number at which 50% of the read lengths are longer than this value", required=True, 
                        type=float) # check
    parser.add_argument("--aligned_read_length_median", help="median read length", required=True, type=float)
    parser.add_argument("--read_qual_median", help="median measure of the uncertainty of base calls", required=True, type=float)

    # Drug Resistance
    parser.add_argument("--drug_resistance_text", help="path of text file used for determining and displaying drug resistances", default=None)
    parser.add_argument("--HRP2", help="value denoting whether the HRP2 marker is present or not -- true or false", default="N/A")
    parser.add_argument("--HRP3", help="value denoting whether the HRP3 marker is present or not -- true or false", default="N/A")

    # Map
    parser.add_argument("--longitude", help="longitude value of where the sample was collected", type=float, default=0)
    parser.add_argument("--latitude", help="latitude value of where the sample collected", type=float, default=0)
    parser.add_argument("--location", help="location where the sample was collected", default="Unknown")
    
    # QC Status
    parser.add_argument("--qc_status", help="status to determine whether or not the sequencing run passes quality control standards", required=True)

    ''' Analysis Page '''
    # required inputs

    # Active Channels
    parser.add_argument("--active_channels", help="number of channels active in the sequencing device", required=True)

    # Q-Scores Plot
    parser.add_argument("--num_reads_q5", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 3", required=True)
    parser.add_argument("--num_reads_q7", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 5", required=True)
    parser.add_argument("--num_reads_q10", help="the number of reads where the probability of a given base call being wrong is 1 in 10", required=True)
    parser.add_argument("--num_reads_q12", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 16", required=True)
    parser.add_argument("--num_reads_q15", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 32", required=True)

    # Sequencing Summary
    parser.add_argument("--sample_prep", help="type of preparation used for the sample", default="N/A")
    parser.add_argument("--analysis_success", help="whether the analysis process completed successfully", required=True)
    parser.add_argument("--aligned_bases", help="total number of bases aligned to the reference genome", required=True)
    parser.add_argument("--aligned_reads", help="total number of reads aligned to the reference genome", required=True)
    parser.add_argument("--fraction_aligned_bases", help="number of bases aligned out of all bases sequenced", required=True, type=float)
    parser.add_argument("--average_identity", help="", required=True, type=float) # check
    
    # Coverage Plot -- in progress
    #parser.add_argument("--coverage_dir", help="location of bed files used for coverage plot")
    parser.add_argument("--fastqc_path", help="location of fastqc_report file; used to locate BAM files for coverage report generation")
    parser.add_argument("--coverage_bin_size", help="number to use as size of bins for coverage plot generation; default is 1500", type=int)

    # parse given arguments
    args = parser.parse_args()
    arg_dict = vars(args)





    '''* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    Prepare arguments for report generation
    '''
    
    # first : summary page
    sample_name = arg_dict['sample_name']

    upload_date = arg_dict['upload_date'][0]
    species = ' '.join(arg_dict['species'])

    info = [upload_date, species, arg_dict['aligned_coverage'], arg_dict['aligned_read_length_n50'], 
            arg_dict['aligned_read_length_median'], arg_dict['read_qual_median']]

    # Check if drug resistance report is provided
    if not arg_dict['drug_resistance_text'] or arg_dict['drug_resistance_text'] == "None":
        resistances = create_drug_table(None)
    else:
        resistances = create_drug_table(arg_dict['drug_resistance_text'])

    qc_status = arg_dict['qc_status']

    # Check if location is given
    if not arg_dict['location']:
        location = "Unknown"
    else:
        location = arg_dict['location']

    # Set default values for map and location info
    make_default = True
    location = "Unknown"
    latitude, longitude = 0, 0
        
    # Check if values for map are provided
    #if arg_dict['latitude'] is not None and arg_dict['longitude'] is not None and arg_dict['location'] is not None:
    #    make_default = False
    #    location = arg_dict['location']
    #    latitude = arg_dict['latitude']
    #    longitude = arg_dict['longitude']

    #if arg_dict['location']: 
    #    location = arg_dict['location']

    #location_info = [latitude, longitude, location]
    #coordinates = [latitude, longitude]

    location_info = [arg_dict['latitude'], arg_dict['longitude'], arg_dict['location']]
    coordinates = [arg_dict['latitude'], arg_dict['longitude']]
    _map = create_map(coordinates, sample_name)

    HRP2 = arg_dict['HRP2']
    HRP3 = arg_dict['HRP3']

    # second : analysis page
    sequencing_summary = [arg_dict['sample_prep'], arg_dict['analysis_success'], arg_dict['aligned_bases'], arg_dict['aligned_reads'], 
                          arg_dict['fraction_aligned_bases'], arg_dict['average_identity']]

    active_channels = arg_dict['active_channels']

    qscorex = [5, 7, 10, 12, 15] # available q-score measures are predetermined
    qscorey = [arg_dict['num_reads_q5'], arg_dict['num_reads_q7'], arg_dict['num_reads_q10'], arg_dict['num_reads_q12'], arg_dict['num_reads_q15']]

    # Create coverage plot and convert it to base64
    coverage_bin_size = arg_dict["coverage_bin_size"]
    # coverage_dir = arg_dict['coverage_dir']
    coverage_plot = plot_coverage("/report-files/data/coverage", sample_name, coverage_bin_size) # default bin size = 1500
    coverage_b64 = img_to_b64(coverage_plot)
    
    # For debugging purposes, save plot:
    coverage_plot.savefig("coverage_plot.jpeg")
    
    # Create summary and analysis objects to be passed 
    summary = Sample(sample_name, HRP2, HRP3, qc_status, resistances, info, _map, location_info)

    analysis = Analysis(sequencing_summary, qscorey, qscorex, active_channels, coverage_b64)

    # Finally, call function to populate and generate the report pages
    create_report(summary, analysis)




