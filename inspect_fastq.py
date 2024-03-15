#!/home/bdhingpit/miniconda3/envs/my_bioinfo_env/bin/python

from chart_studio import plotly as py
import plotly.graph_objects as go
import plotly.express as px
import colorlover as cl
import skbio
from skbio.sequence import DNA, RNA
import pandas as pd
import itertools
import numpy as np
import argparse
import time
from warnings import simplefilter

start_time = time.time()

#Module: per_base_qc_plot
#Module: fastq_stat
def define_layout(plot_title, x_title, y_title, y_range, show_legend):
	layout = go.Layout(
		title=dict(
			text=plot_title,
			x=0.5
		),
		template="plotly_white",
		yaxis=dict(
			title=y_title,
			showgrid=True,
			gridcolor="#d9d4d3",
			zerolinecolor="#d9d4d3",
			range=y_range
		),
		xaxis=dict(
			title=x_title,
		),
		font=dict(
			family="Open Sans", 
			size=16, 
			color="#2e1c18"
		),
		showlegend=show_legend
	)

	return layout

#Module: per_base_qc_plot
def module_per_base_qc_plot(args):
	fastq_df = tabulate_fastq_info(args.fastq, args.encoding, args.num_seqs)
	show_per_base_qc_plot(fastq_df, None)

	return None

#Module: per_base_qc_plot
def tabulate_fastq_info(fastq, encoding, num_seqs):
	imported_fastq = list(skbio.io.read(fastq, format='fastq', verify=False, variant=encoding))

	fastq_df = pd.DataFrame()

	for count, seq in enumerate(imported_fastq):
		fastq_df[count] = seq.positional_metadata.quality

		if count==num_seqs-1:
			break

	fastq_df.index +=1

	return fastq_df

#Module: per_base_qc_plot
def define_color_scale():
	col_scales = cl.scales['4']['div']['RdYlGn']
	col_scales_40 = cl.interp(col_scales, 40)

	return col_scales_40

#Module: per_base_qc_plot
def create_boxes(fastq_df, col_scales):
	traces = []

	for base in range(len(fastq_df)):
		traces.append(
			go.Box(
				#y=fastq_df.iloc[base].values,
				name="Base Position Quality",
				x=[base+1],
				boxpoints=False,
				whiskerwidth=0.5,
				marker=dict(
					size=.1,
					color=col_scales[int(round(fastq_df.iloc[base].mean(), 0))]
		        ),
		        line=dict(width=1),
				q1=[np.nanpercentile(fastq_df.iloc[base].values, 25)],
				q3=[np.nanpercentile(fastq_df.iloc[base].values, 75)],
        		median=[np.nanpercentile(fastq_df.iloc[base].values, 50)],
        		lowerfence=[np.nanpercentile(fastq_df.iloc[base].values, 10)],
        		upperfence=[np.nanpercentile(fastq_df.iloc[base].values, 90)],
        		hoverlabel=dict(
        			namelength=-1,
        			align="left"
        		)
	        )
		)

	return traces

#Module: per_base_qc_plot
def show_per_base_qc_plot(fastq_df, output):
	col_scales_40 = define_color_scale()
	traces = create_boxes(fastq_df, col_scales_40)
	layout = define_layout(
		"Per-Base Quality Score", 
		"Base Position", 
		"Quality Score",
		[0, 41], #Max possible Q-score = 41
		False
	)

	fig = go.Figure(data=traces, layout=layout)
	#fig.update_yaxes(ticksuffix = "    ")

	fig.show()

#Module: fastq_stat
def module_fastq_stat(args):
	imported_fastq = list(skbio.io.read(args.fastq, format='fastq', verify=False, variant=args.encoding))

	show_fastq_stat_plots(imported_fastq)
	#get_gc_content(imported_fastq)

	return None

#Module: fastq_stat
def show_fastq_stat_plots(fastq):
	hist_trace = create_seq_length_hist(fastq)

	fig = go.Figure(data=traces)

	fig.show()

#Module: fastq_stat
def create_seq_length_hist(fastq):
	fastq_seq_lengths = [len(seq) for seq in fastq]

	hist_trace = go.Histogram(
		x=fastq_seq_lengths
	)

	return hist_trace

#Module: fastq_stat
def create_gc_content(fastq):
	concat_seq = ""

	for seq in fastq:
		concat_seq += str(seq)

	gc_content = DNA(concat_seq).gc_content()

#Module: fastq_stat
def get_num_seqs():
	len(imported_fastq)




def main(args):
	#Ignore pandas warnings
	simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
	
	args.func(args)



if __name__ == "__main__":
	parent_parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		description="Inspect FASTQ files interactively",
		add_help=False
	)

	#Define arguments always required
	required_args = parent_parser.add_argument_group('Required arguments')

	required_args.add_argument(
		"--fastq", 
		dest="fastq", 
		type=str,
		required=True,
		metavar="STRING",
		help="Path to FASTQ file"
	)

	required_args.add_argument(
		"--encoding", 
		dest="encoding", 
		type=str,
		required=False,
		choices=["sanger", "illumina1.3", "illumina1.8", "solexa"],
		default="illumina1.8",
		metavar="STRING",
		help="Encoding system used in FASTQ file. Default: \"illumina1.8\" || Choices: [\"sanger\", \"illumina1.3\", \"illumina1.8\", \"solexa\"]"
	)

	#Define main parser
	main_parser = argparse.ArgumentParser()

	subparser = main_parser.add_subparsers(
		title="Analysis type",
		dest="analysis_type",
		description="Available subcommands"
	)

	#Parser for per_base_qc_plot module
	per_base_qc_plot_subparser = subparser.add_parser("per_base_qc_plot", parents=[required_args])

	per_base_qc_plot_subparser.add_argument(
		"--num_seqs", 
		dest="num_seqs", 
		type=int,
		required=False,
		default=10000,
		metavar="INTEGER",
		help="Number of sequences to subsample for per-base QC plot. Default: 10000"
	)

	per_base_qc_plot_subparser.set_defaults(func=module_per_base_qc_plot)

	#Parser for fastq_stat module
	fastq_stat_subparser = subparser.add_parser("fastq_stat", parents=[required_args])

	fastq_stat_subparser.set_defaults(func=module_fastq_stat)


	args = main_parser.parse_args()

	main(args)