#!/usr/bin/env python3
"""
BELT Mutation Analysis Visualizer
Standalone Plotly Dash web application for visualizing mutation analysis results
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import logging

# Plotly imports
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, callback, dash_table
import dash_bootstrap_components as dbc

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "BELT Mutation Analysis Visualizer"

# Layout will be defined after data loading


# Global variables to store data
analysis_data = None
summary_df = None

# Auto-load data on startup
def load_data_on_startup():
    """Load data automatically on startup"""
    global analysis_data, summary_df
    
    data_dir = "/Users/anitido/Documents/GitHub/belt/belt/Belt Pipeline/data/9-mutation_analysis"
    try:
        data_path = Path(data_dir)
        if data_path.exists():
            # Load the main analysis report
            report_file = data_path / "mutation_analysis_report.json"
            if report_file.exists():
                with open(report_file, 'r') as f:
                    analysis_data = json.load(f)
                
                # Load CSV files if they exist
                summary_file = data_path / "mutation_summary.csv"
                if summary_file.exists():
                    summary_df = pd.read_csv(summary_file)
                
                logger.info(f"Auto-loaded data: {len(analysis_data['results'])} samples found.")
                logger.info(f"Sample names: {[r['sample_name'] for r in analysis_data['results']]}")
                return True
    except Exception as e:
        logger.error(f"Error auto-loading data: {e}")
    
    return False

# Load data on startup
data_loaded = load_data_on_startup()

# Create a simple initial content function
def create_simple_summary():
    """Create a simple summary when data is loaded"""
    if not analysis_data:
        return html.Div([
            html.H4("No Data Available"),
            html.P("Please ensure the mutation analysis data exists in the specified directory.")
        ])
    
    # Create basic summary
    total_samples = len(analysis_data['results'])
    total_consensus = sum(r['consensus_count'] for r in analysis_data['results'])
    total_mutations = sum(r['summary_stats']['total_mutations'] for r in analysis_data['results'])
    
    # Create sortable sample details table
    sample_data = []
    for r in analysis_data['results']:
        sample_data.append({
            'Sample': r['sample_name'],
            'Consensus Count': r['consensus_count'],
            'Total Mutations': r['summary_stats']['total_mutations'],
            'Avg Mutation Rate': r['summary_stats']['average_mutation_rate'],
            'Substitutions': r['summary_stats']['mutation_types']['substitution'],
            'Insertions': r['summary_stats']['mutation_types']['insertion'],
            'Deletions': r['summary_stats']['mutation_types']['deletion']
        })
    
    sample_df = pd.DataFrame(sample_data)
    
    return html.Div([
        html.H4("📊 Analysis Summary"),
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H4(f"{total_samples}", className="text-primary"),
                        html.P("Samples", className="mb-0")
                    ])
                ])
            ], width=3),
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H4(f"{total_consensus}", className="text-success"),
                        html.P("Consensus Sequences", className="mb-0")
                    ])
                ])
            ], width=3),
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H4(f"{total_mutations}", className="text-warning"),
                        html.P("Total Mutations", className="mb-0")
                    ])
                ])
            ], width=3),
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H4(f"{total_mutations/total_consensus:.2f}", className="text-info"),
                        html.P("Avg Mutations/Read", className="mb-0")
                    ])
                ])
            ], width=3)
        ], className="mb-4"),
        
        html.H5("Sample Details (Sortable)"),
        dash_table.DataTable(
            data=sample_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Consensus Count", "id": "Consensus Count", "type": "numeric"},
                {"name": "Total Mutations", "id": "Total Mutations", "type": "numeric"},
                {"name": "Avg Mutation Rate", "id": "Avg Mutation Rate", "type": "numeric", "format": {"specifier": ".6f"}},
                {"name": "Substitutions", "id": "Substitutions", "type": "numeric"},
                {"name": "Insertions", "id": "Insertions", "type": "numeric"},
                {"name": "Deletions", "id": "Deletions", "type": "numeric"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
        )
    ])


def create_rates_tab():
    """Create mutation rates visualization tab"""
    if not analysis_data:
        return html.Div("Please load data first.")
    
    # Create mutation rate visualization
    samples = [r['sample_name'] for r in analysis_data['results']]
    avg_rates = [r['summary_stats']['average_mutation_rate'] for r in analysis_data['results']]
    median_rates = [r['summary_stats']['median_mutation_rate'] for r in analysis_data['results']]
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name='Average Rate',
        x=samples,
        y=avg_rates,
        marker_color='lightblue'
    ))
    fig.add_trace(go.Bar(
        name='Median Rate',
        x=samples,
        y=median_rates,
        marker_color='lightcoral'
    ))
    
    fig.update_layout(
        title="Mutation Rates by Sample",
        xaxis_title="Sample",
        yaxis_title="Mutation Rate",
        barmode='group'
    )
    
    # Create sortable rate statistics table
    rate_data = []
    for r in analysis_data['results']:
        rate_data.append({
            'Sample': r['sample_name'],
            'Average Rate': r['summary_stats']['average_mutation_rate'],
            'Median Rate': r['summary_stats']['median_mutation_rate'],
            'Std Dev': r['summary_stats']['mutation_rate_std']
        })
    
    rate_df = pd.DataFrame(rate_data)
    
    return html.Div([
        html.H4("📈 Mutation Rates Analysis"),
        dcc.Graph(figure=fig),
        html.H5("Rate Statistics (Sortable)"),
        dash_table.DataTable(
            data=rate_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Average Rate", "id": "Average Rate", "type": "numeric", "format": {"specifier": ".6f"}},
                {"name": "Median Rate", "id": "Median Rate", "type": "numeric", "format": {"specifier": ".6f"}},
                {"name": "Std Dev", "id": "Std Dev", "type": "numeric", "format": {"specifier": ".6f"}}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ]
        )
    ])


def create_per_read_tab():
    """Create per-read mutations analysis tab"""
    if not analysis_data:
        return html.Div("Please load data first.")
    
    # Create mutations per read visualization
    all_mutations_per_read = []
    for result in analysis_data['results']:
        for read_data in result['per_read_analysis']:
            all_mutations_per_read.append(read_data['mutations_count'])
    
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=all_mutations_per_read,
        nbinsx=20,
        marker_color='lightgreen'
    ))
    
    fig.update_layout(
        title="Distribution of Mutations per Read",
        xaxis_title="Number of Mutations",
        yaxis_title="Frequency"
    )
    
    # Create sortable per-read statistics table
    per_read_data = []
    for r in analysis_data['results']:
        per_read_data.append({
            'Sample': r['sample_name'],
            'Mean Mutations/Read': r['summary_stats']['mutations_per_read_mean'],
            'Std Dev': r['summary_stats']['mutations_per_read_std'],
            'Total Reads': r['consensus_count']
        })
    
    per_read_df = pd.DataFrame(per_read_data)
    
    # Create detailed per-read data table
    detailed_per_read_data = []
    for result in analysis_data['results']:
        sample_name = result['sample_name']
        for read_data in result['per_read_analysis']:
            detailed_per_read_data.append({
                'Sample': sample_name,
                'Consensus ID': read_data['consensus_id'],
                'Mutations Count': read_data['mutations_count'],
                'Consensus Length': read_data['consensus_length'],
                'Alignment Score': read_data['alignment_score'],
                'Matches': read_data['matches'],
                'Mismatches': read_data['mismatches'],
                'Gaps': read_data['gaps']
            })
    
    detailed_df = pd.DataFrame(detailed_per_read_data)
    
    return html.Div([
        html.H4("📋 Per-Read Mutations Analysis"),
        dcc.Graph(figure=fig),
        html.H5("Per-Read Statistics (Sortable)"),
        dash_table.DataTable(
            data=per_read_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Mean Mutations/Read", "id": "Mean Mutations/Read", "type": "numeric", "format": {"specifier": ".2f"}},
                {"name": "Std Dev", "id": "Std Dev", "type": "numeric", "format": {"specifier": ".2f"}},
                {"name": "Total Reads", "id": "Total Reads", "type": "numeric"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
        ),
        html.H5("Detailed Per-Read Data (Sortable)"),
        dash_table.DataTable(
            data=detailed_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Consensus ID", "id": "Consensus ID", "type": "text"},
                {"name": "Mutations Count", "id": "Mutations Count", "type": "numeric"},
                {"name": "Consensus Length", "id": "Consensus Length", "type": "numeric"},
                {"name": "Alignment Score", "id": "Alignment Score", "type": "numeric", "format": {"specifier": ".1f"}},
                {"name": "Matches", "id": "Matches", "type": "numeric"},
                {"name": "Mismatches", "id": "Mismatches", "type": "numeric"},
                {"name": "Gaps", "id": "Gaps", "type": "numeric"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ]
        )
    ])


def create_distribution_tab():
    """Create mutation type distribution tab"""
    if not analysis_data:
        return html.Div("Please load data first.")
    
    # Get reference sequence to calculate nucleotide composition
    reference_path = "/Users/anitido/Documents/GitHub/belt/belt/Belt Pipeline/data/0-reference_sequence"
    reference_sequence = None
    nucleotide_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    
    try:
        for result in analysis_data['results']:
            if 'reference_id' in result:
                ref_id = result['reference_id']
                ref_file = f"{reference_path}/{ref_id}.fa"
                if os.path.exists(ref_file):
                    with open(ref_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 2:
                            reference_sequence = lines[1].strip().upper()
                            # Count nucleotides in reference
                            for base in reference_sequence:
                                if base in nucleotide_counts:
                                    nucleotide_counts[base] += 1
                            break
    except Exception as e:
        logger.warning(f"Could not load reference sequence: {e}")
    
    # Collect all mutations and create mutation type distribution
    all_mutations = []
    mutation_type_counts = {}
    mutation_type_rates = {}
    
    for result in analysis_data['results']:
        sample_name = result['sample_name']
        for mutation in result['summary_stats']['mutations_by_position']:
            ref_base = mutation['ref_base']
            cons_base = mutation['cons_base']
            mutation_type = f"{ref_base}→{cons_base}"
            
            all_mutations.append({
                'Sample': sample_name,
                'Position': mutation['position'],
                'Type': mutation['type'],
                'Ref Base': ref_base,
                'Consensus Base': cons_base,
                'Mutation Type': mutation_type
            })
            
            # Count mutation types
            if mutation_type not in mutation_type_counts:
                mutation_type_counts[mutation_type] = 0
            mutation_type_counts[mutation_type] += 1
    
    # Calculate rates based on reference nucleotide composition
    total_reference_bases = sum(nucleotide_counts.values())
    for mutation_type, count in mutation_type_counts.items():
        ref_base = mutation_type.split('→')[0]
        if ref_base in nucleotide_counts and nucleotide_counts[ref_base] > 0:
            mutation_type_rates[mutation_type] = count / nucleotide_counts[ref_base]
        else:
            mutation_type_rates[mutation_type] = 0
    
    # Create mutation type distribution visualization
    mutation_types = list(mutation_type_counts.keys())
    counts = list(mutation_type_counts.values())
    rates = [mutation_type_rates[mt] for mt in mutation_types]
    
    # Create bar chart for mutation type counts
    fig_counts = go.Figure()
    fig_counts.add_trace(go.Bar(
        x=mutation_types,
        y=counts,
        marker_color='lightblue',
        name='Counts'
    ))
    
    fig_counts.update_layout(
        title="Mutation Type Distribution - Total Counts",
        xaxis_title="Mutation Type",
        yaxis_title="Count",
        showlegend=False
    )
    
    # Create bar chart for mutation type rates
    fig_rates = go.Figure()
    fig_rates.add_trace(go.Bar(
        x=mutation_types,
        y=rates,
        marker_color='lightcoral',
        name='Rates'
    ))
    
    fig_rates.update_layout(
        title="Mutation Type Distribution - Rates (per reference nucleotide)",
        xaxis_title="Mutation Type",
        yaxis_title="Rate",
        showlegend=False
    )
    
    # Create pie chart for mutation type distribution
    fig_pie = go.Figure()
    fig_pie.add_trace(go.Pie(
        labels=mutation_types,
        values=counts,
        name="Mutation Types"
    ))
    
    fig_pie.update_layout(
        title="Mutation Type Distribution - Pie Chart"
    )
    
    # Create mutation type summary table
    mutation_summary_data = []
    for mutation_type in mutation_types:
        mutation_summary_data.append({
            'Mutation Type': mutation_type,
            'Count': mutation_type_counts[mutation_type],
            'Rate': mutation_type_rates[mutation_type],
            'Reference Count': nucleotide_counts[mutation_type.split('→')[0]] if mutation_type.split('→')[0] in nucleotide_counts else 0
        })
    
    mutation_summary_df = pd.DataFrame(mutation_summary_data)
    mutation_summary_df = mutation_summary_df.sort_values('Count', ascending=False)
    
    # Create detailed mutation table
    mutation_df = pd.DataFrame(all_mutations)
    
    return html.Div([
        html.H4("🗺️ Mutation Type Distribution Analysis"),
        html.P(f"Total reference length: {total_reference_bases} bases"),
        html.P(f"Reference composition: A={nucleotide_counts['A']}, T={nucleotide_counts['T']}, G={nucleotide_counts['G']}, C={nucleotide_counts['C']}"),
        
        # Mutation type distribution charts
        dbc.Row([
            dbc.Col([
                dcc.Graph(figure=fig_counts)
            ], width=6),
            dbc.Col([
                dcc.Graph(figure=fig_rates)
            ], width=6)
        ], className="mb-4"),
        
        dbc.Row([
            dbc.Col([
                dcc.Graph(figure=fig_pie)
            ], width=12)
        ], className="mb-4"),
        
        html.H5("Mutation Type Summary (Sortable)"),
        dash_table.DataTable(
            data=mutation_summary_df.to_dict('records'),
            columns=[
                {"name": "Mutation Type", "id": "Mutation Type", "type": "text"},
                {"name": "Count", "id": "Count", "type": "numeric"},
                {"name": "Rate", "id": "Rate", "type": "numeric", "format": {"specifier": ".6f"}},
                {"name": "Reference Count", "id": "Reference Count", "type": "numeric"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
        ),
        
        html.H5("Detailed Mutation Data (Sortable)"),
        dash_table.DataTable(
            data=mutation_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Position", "id": "Position", "type": "numeric"},
                {"name": "Type", "id": "Type", "type": "text"},
                {"name": "Ref Base", "id": "Ref Base", "type": "text"},
                {"name": "Consensus Base", "id": "Consensus Base", "type": "text"},
                {"name": "Mutation Type", "id": "Mutation Type", "type": "text"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ]
        )
    ])


def create_comparison_tab():
    """Create sample comparison tab"""
    if not analysis_data:
        return html.Div("Please load data first.")
    
    # Create comparison visualization
    samples = [r['sample_name'] for r in analysis_data['results']]
    mutations = [r['summary_stats']['total_mutations'] for r in analysis_data['results']]
    consensus = [r['consensus_count'] for r in analysis_data['results']]
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=consensus,
        y=mutations,
        mode='markers+text',
        text=samples,
        textposition="top center",
        marker=dict(size=15, color='red')
    ))
    
    fig.update_layout(
        title="Sample Comparison: Consensus vs Mutations",
        xaxis_title="Number of Consensus Sequences",
        yaxis_title="Total Mutations"
    )
    
    # Create sortable comparison summary table
    comparison_data = []
    for r in analysis_data['results']:
        comparison_data.append({
            'Sample': r['sample_name'],
            'Consensus Count': r['consensus_count'],
            'Total Mutations': r['summary_stats']['total_mutations'],
            'Mutation Rate': r['summary_stats']['average_mutation_rate'],
            'Substitutions': r['summary_stats']['mutation_types']['substitution'],
            'Insertions': r['summary_stats']['mutation_types']['insertion'],
            'Deletions': r['summary_stats']['mutation_types']['deletion']
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    return html.Div([
        html.H4("🔍 Sample Comparison Analysis"),
        dcc.Graph(figure=fig),
        html.H5("Comparison Summary (Sortable)"),
        dash_table.DataTable(
            data=comparison_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "Sample", "type": "text"},
                {"name": "Consensus Count", "id": "Consensus Count", "type": "numeric"},
                {"name": "Total Mutations", "id": "Total Mutations", "type": "numeric"},
                {"name": "Mutation Rate", "id": "Mutation Rate", "type": "numeric", "format": {"specifier": ".6f"}},
                {"name": "Substitutions", "id": "Substitutions", "type": "numeric"},
                {"name": "Insertions", "id": "Insertions", "type": "numeric"},
                {"name": "Deletions", "id": "Deletions", "type": "numeric"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ]
        )
    ])


def create_raw_data_tab():
    """Create raw data viewer tab"""
    if not analysis_data:
        return html.Div("Please load data first.")
    
    return html.Div([
        html.H4("📄 Raw Data Viewer"),
        html.H5("JSON Data Structure"),
        html.Pre(
            json.dumps(analysis_data, indent=2),
            style={
                'background-color': '#f8f9fa',
                'border': '1px solid #dee2e6',
                'border-radius': '0.375rem',
                'padding': '1rem',
                'max-height': '500px',
                'overflow-y': 'auto'
            }
        )
    ])

# Create initial content based on data loading status
if data_loaded and analysis_data:
    initial_status = dbc.Alert(f"✅ Data auto-loaded successfully! {len(analysis_data['results'])} samples found.", color="success")
    initial_content = create_simple_summary()
else:
    initial_status = dbc.Alert("❌ No data found. Please check the data directory.", color="warning")
    initial_content = html.Div([
        html.H4("No Data Available"),
        html.P("Please ensure the mutation analysis data exists in the specified directory.")
    ])

# Update the layout to show initial content
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("🧬 BELT Mutation Analysis Visualizer", className="text-center mb-4"),
            html.Hr(),
            
            # Data loading section
            dbc.Card([
                dbc.CardBody([
                    html.H4("📊 Data Loading", className="card-title"),
                    html.P("Load mutation analysis results from your BELT pipeline output."),
                    dbc.Row([
                        dbc.Col([
                            html.Label("Data Directory:"),
                            dcc.Input(
                                id="data-directory",
                                type="text",
                                placeholder="/path/to/mutation_analysis/",
                                value="/Users/anitido/Documents/GitHub/belt/belt/Belt Pipeline/data/9-mutation_analysis",
                                style={"width": "100%"}
                            )
                        ], width=8),
                        dbc.Col([
                            html.Br(),
                            dbc.Button("Load Data", id="load-data-btn", color="primary", className="mt-2")
                        ], width=4)
                    ]),
                    html.Div(id="data-status", className="mt-2"),
                    html.Div(id="auto-status", className="mt-2", children=initial_status)
                ])
            ], className="mb-4"),
            
            # Navigation tabs
            dbc.Tabs([
                dbc.Tab(label="📊 Summary Overview", tab_id="summary", id="summary-tab"),
                dbc.Tab(label="📈 Mutation Rates", tab_id="rates", id="rates-tab"),
                dbc.Tab(label="📋 Per-Read Analysis", tab_id="per-read", id="per-read-tab"),
                dbc.Tab(label="🗺️ Mutation Distribution", tab_id="distribution", id="distribution-tab"),
                dbc.Tab(label="🔍 Sample Comparison", tab_id="comparison", id="comparison-tab"),
                dbc.Tab(label="🧬 Alignment Visualizer", tab_id="alignment", id="alignment-tab"),
                dbc.Tab(label="📄 Raw Data", tab_id="raw-data", id="raw-data-tab")
            ], id="main-tabs", active_tab="summary"),
            
            # Tab content
            html.Div(id="tab-content", children=initial_content)
        ])
    ])
], fluid=True)

@callback(
    [Output("data-status", "children"),
     Output("tab-content", "children")],
    [Input("load-data-btn", "n_clicks"),
     Input("main-tabs", "active_tab")],
    [dash.dependencies.State("data-directory", "value")]
)
def update_content(load_clicks, active_tab, data_dir):
    """Handle both data loading and tab switching"""
    global analysis_data, summary_df
    
    # Check if this is a data loading trigger
    if load_clicks and data_dir and data_dir.strip():
        try:
            data_path = Path(data_dir)
            if not data_path.exists():
                return dbc.Alert(f"Directory does not exist: {data_dir}", color="danger"), ""
            
            # Load the main analysis report
            report_file = data_path / "mutation_analysis_report.json"
            if not report_file.exists():
                return dbc.Alert(f"mutation_analysis_report.json not found in directory: {data_dir}", color="danger"), ""
            
            with open(report_file, 'r') as f:
                analysis_data = json.load(f)
            
            # Load CSV files if they exist
            summary_file = data_path / "mutation_summary.csv"
            if summary_file.exists():
                summary_df = pd.read_csv(summary_file)
            else:
                summary_df = None
            
            # Debug info
            logger.info(f"Data loaded successfully! {len(analysis_data['results'])} samples found.")
            logger.info(f"Analysis data keys: {list(analysis_data.keys())}")
            logger.info(f"Sample names: {[r['sample_name'] for r in analysis_data['results']]}")
            
            status = dbc.Alert(f"✅ Data loaded successfully! {len(analysis_data['results'])} samples found.", 
                              color="success")
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            return dbc.Alert(f"❌ Error loading data: {str(e)}", color="danger"), ""
    else:
        status = ""
    
    # Handle tab content
    if not analysis_data:
        return status, html.Div("Please load data first.")
    
    if active_tab == "summary":
        content = create_simple_summary()
    elif active_tab == "rates":
        content = create_rates_tab()
    elif active_tab == "per-read":
        content = create_per_read_tab()
    elif active_tab == "distribution":
        content = create_distribution_tab()
    elif active_tab == "comparison":
        content = create_comparison_tab()
    elif active_tab == "alignment":
        content = create_alignment_visualizer_tab()
    elif active_tab == "raw-data":
        content = create_raw_data_tab()
    else:
        content = create_simple_summary()
    
    return status, content


def create_summary_tab():
    """Create the summary overview tab content"""
    global analysis_data, summary_df
    
    if not analysis_data:
        return html.Div("Please load data first.")
    
    data = analysis_data
    results = data['results']
    
    # Calculate overall statistics
    total_samples = len([r for r in results if "error" not in r])
    total_mutations = sum(r.get("summary_stats", {}).get("total_mutations", 0) for r in results if "error" not in r)
    total_consensus = sum(r.get("consensus_count", 0) for r in results if "error" not in r)
    
    # Create summary cards
    cards = [
        dbc.Card([
            dbc.CardBody([
                html.H4(f"{total_samples}", className="card-title"),
                html.P("Samples Analyzed", className="card-text")
            ])
        ], color="primary", outline=True),
        dbc.Card([
            dbc.CardBody([
                html.H4(f"{total_consensus:,}", className="card-title"),
                html.P("Consensus Sequences", className="card-text")
            ])
        ], color="info", outline=True),
        dbc.Card([
            dbc.CardBody([
                html.H4(f"{total_mutations:,}", className="card-title"),
                html.P("Total Mutations", className="card-text")
            ])
        ], color="warning", outline=True)
    ]
    
    # Sample table
    if summary_df is not None:
        sample_table = dash_table.DataTable(
            data=summary_df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in summary_df.columns],
            style_cell={'textAlign': 'left'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
            page_size=10
        )
    else:
        sample_table = html.P("Summary data not available")
    
    return html.Div([
        dbc.Row([dbc.Col(card, width=4) for card in cards], className="mb-4"),
        html.H4("Sample Summary"),
        sample_table
    ])


def create_mutation_rates_tab():
    """Create the mutation rates visualization tab"""
    global summary_df
    
    if summary_df is None:
        return html.Div("Summary data not available. Please ensure mutation_summary.csv exists.")
    
    df = summary_df
    
    # Create mutation rate comparison chart
    fig_rates = px.bar(
        df, 
        x='sample_name', 
        y='average_mutation_rate',
        title="Average Mutation Rate by Sample",
        labels={'average_mutation_rate': 'Mutation Rate', 'sample_name': 'Sample'}
    )
    fig_rates.update_layout(showlegend=False)
    
    # Create mutation rate distribution
    fig_dist = px.histogram(
        df, 
        x='average_mutation_rate',
        title="Distribution of Mutation Rates",
        labels={'average_mutation_rate': 'Mutation Rate', 'count': 'Number of Samples'}
    )
    
    return html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Graph(figure=fig_rates)
            ], width=6),
            dbc.Col([
                dcc.Graph(figure=fig_dist)
            ], width=6)
        ])
    ])


def create_per_read_tab():
    """Create the mutations per read visualization tab"""
    global analysis_data
    
    if not analysis_data:
        return html.Div("Analysis data not available.")
    
    # Collect per-read data from all samples
    per_read_data = []
    for result in analysis_data['results']:
        if "error" not in result and "per_read_analysis" in result:
            sample_name = result["sample_name"]
            for read in result["per_read_analysis"]:
                per_read_data.append({
                    "sample": sample_name,
                    "consensus_id": read["consensus_id"],
                    "mutation_count": read["mutation_count"],
                    "consensus_length": read["consensus_length"],
                    "mutation_rate": read["mutation_count"] / read["consensus_length"] if read["consensus_length"] > 0 else 0
                })
    
    if not per_read_data:
        return html.Div("No per-read data available.")
    
    df = pd.DataFrame(per_read_data)
    
    # Create mutations per read histogram
    fig_hist = px.histogram(
        df, 
        x='mutation_count',
        title="Distribution of Mutations per Read",
        labels={'mutation_count': 'Mutations per Read', 'count': 'Number of Reads'}
    )
    
    # Create box plot by sample
    fig_box = px.box(
        df, 
        x='sample', 
        y='mutation_count',
        title="Mutations per Read by Sample",
        labels={'mutation_count': 'Mutations per Read', 'sample': 'Sample'}
    )
    
    return html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Graph(figure=fig_hist)
            ], width=6),
            dbc.Col([
                dcc.Graph(figure=fig_box)
            ], width=6)
        ])
    ])


def create_comparison_tab():
    """Create the sample comparison tab"""
    global summary_df
    
    if summary_df is None:
        return html.Div("Summary data not available.")
    
    df = summary_df
    
    # Create comparison scatter plot
    fig_scatter = px.scatter(
        df,
        x='total_mutations',
        y='average_mutation_rate',
        size='consensus_count',
        color='sample_name',
        title="Sample Comparison: Total Mutations vs Mutation Rate",
        labels={'total_mutations': 'Total Mutations', 'average_mutation_rate': 'Average Mutation Rate'}
    )
    
    # Create mutation type comparison
    mutation_types = ['substitutions', 'insertions', 'deletions']
    fig_types = go.Figure()
    
    for i, mut_type in enumerate(mutation_types):
        fig_types.add_trace(go.Bar(
            name=mut_type.title(),
            x=df['sample_name'],
            y=df[mut_type],
            offsetgroup=i
        ))
    
    fig_types.update_layout(
        title="Mutation Types by Sample",
        xaxis_title="Sample",
        yaxis_title="Count"
    )
    
    return html.Div([
        dbc.Row([
            dbc.Col([
                dcc.Graph(figure=fig_scatter)
            ], width=6),
            dbc.Col([
                dcc.Graph(figure=fig_types)
            ], width=6)
        ])
    ])


def create_alignment_visualizer_tab():
    """Create alignment visualizer tab showing consensus reads aligned to reference with mutations highlighted"""
    global analysis_data
    
    if not analysis_data:
        return html.Div("Analysis data not available.")
    
    # Get reference sequence
    reference_path = "/Users/anitido/Documents/GitHub/belt/belt/Belt Pipeline/data/0-reference_sequence"
    reference_sequence = None
    
    try:
        # Try to find the reference sequence file
        for result in analysis_data['results']:
            if 'reference_id' in result:
                ref_id = result['reference_id']
                ref_file = f"{reference_path}/{ref_id}.fa"
                if os.path.exists(ref_file):
                    with open(ref_file, 'r') as f:
                        lines = f.readlines()
                        if len(lines) >= 2:
                            reference_sequence = lines[1].strip()
                            break
    except Exception as e:
        logger.warning(f"Could not load reference sequence: {e}")
    
    if not reference_sequence:
        return html.Div([
            html.H4("🧬 Alignment Visualizer"),
            html.P("Reference sequence not found. Please ensure the reference file exists.")
        ])
    
    # Create alignment visualization
    alignment_data = []
    for result in analysis_data['results']:
        sample_name = result['sample_name']
        for read_data in result['per_read_analysis']:
            consensus_id = read_data['consensus_id']
            mutations = read_data.get('mutations', [])
            
            # Create mutation positions set for quick lookup
            mutation_positions = {mut['position']: mut for mut in mutations}
            
            # Create alignment row
            alignment_row = {
                'sample': sample_name,
                'consensus_id': consensus_id,
                'mutations_count': len(mutations),
                'alignment_score': read_data.get('alignment_score', 0),
                'matches': read_data.get('matches', 0),
                'mismatches': read_data.get('mismatches', 0),
                'gaps': read_data.get('gaps', 0)
            }
            
            # Add mutation details
            mutation_details = []
            for pos, mut in mutation_positions.items():
                mutation_details.append(f"Pos {pos}: {mut['ref_base']}→{mut['cons_base']} ({mut['type']})")
            
            alignment_row['mutation_details'] = '; '.join(mutation_details) if mutation_details else 'No mutations'
            alignment_data.append(alignment_row)
    
    alignment_df = pd.DataFrame(alignment_data)
    
    # Create interactive alignment plot
    fig = go.Figure()
    
    # Add reference sequence as a horizontal line
    ref_length = len(reference_sequence)
    fig.add_trace(go.Scatter(
        x=[0, ref_length],
        y=[0, 0],
        mode='lines',
        line=dict(color='black', width=4),
        name='Reference Sequence',
        showlegend=True,
        hovertemplate='Reference: %{x}<br>Length: %{x}<extra></extra>'
    ))
    
    # Define mutation type symbols and colors - organized by reference nucleotide
    # Same shape for all mutations from same reference nucleotide
    # Consistent colors for each target nucleotide across all reference nucleotides
    mutation_symbols = {
        # A→X mutations (all circles, consistent target colors)
        'A→G': {'symbol': 'circle', 'color': 'red'},
        'A→T': {'symbol': 'circle', 'color': 'blue'},
        'A→C': {'symbol': 'circle', 'color': 'green'},
        # G→X mutations (all squares, consistent target colors)
        'G→A': {'symbol': 'square', 'color': 'orange'},
        'G→T': {'symbol': 'square', 'color': 'blue'},
        'G→C': {'symbol': 'square', 'color': 'green'},
        # T→X mutations (all diamonds, consistent target colors)
        'T→A': {'symbol': 'diamond', 'color': 'orange'},
        'T→G': {'symbol': 'diamond', 'color': 'red'},
        'T→C': {'symbol': 'diamond', 'color': 'green'},
        # C→X mutations (all triangles, consistent target colors)
        'C→A': {'symbol': 'triangle-up', 'color': 'orange'},
        'C→G': {'symbol': 'triangle-up', 'color': 'red'},
        'C→T': {'symbol': 'triangle-up', 'color': 'blue'}
    }
    
    # Create legend entries for all possible mutation types (hidden traces)
    for mutation_type, symbol_info in mutation_symbols.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None],  # Use None instead of empty arrays
            mode='markers',
            marker=dict(
                color=symbol_info['color'],
                size=12,
                symbol=symbol_info['symbol'],
                line=dict(width=2, color='black')
            ),
            name=mutation_type,
            showlegend=True,
            visible=True  # Make sure they're visible in legend
        ))
    
    # Add consensus reads with mutations highlighted
    y_offset = 0.1  # Start very close to reference
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    
    for i, (_, row) in enumerate(alignment_df.iterrows()):
        # Get mutations for this read
        mutations = []
        for result in analysis_data['results']:
            for read_data in result['per_read_analysis']:
                if read_data['consensus_id'] == row['consensus_id']:
                    mutations = read_data.get('mutations', [])
                    break
        
        # Create mutation positions
        mutation_positions = [mut['position'] for mut in mutations]
        
        # Add consensus read line
        fig.add_trace(go.Scatter(
            x=[0, ref_length],
            y=[y_offset, y_offset],
            mode='lines',
            line=dict(color=colors[i % len(colors)], width=3),
            name=f"Read {i+1} ({row['mutations_count']} muts)",
            showlegend=False,  # Hide read lines from legend
            hovertemplate=f"Read: {row['consensus_id']}<br>Mutations: {row['mutations_count']}<br>Score: {row['alignment_score']}<extra></extra>"
        ))
        
        # Add mutation markers with different symbols for each mutation type
        if mutation_positions:
            # Group mutations by type for separate traces
            mutations_by_type = {}
            for mut in mutations:
                mutation_type = f"{mut['ref_base']}→{mut['cons_base']}"
                if mutation_type not in mutations_by_type:
                    mutations_by_type[mutation_type] = []
                mutations_by_type[mutation_type].append(mut)
            
            # Create separate traces for each mutation type (hidden from legend)
            for mutation_type, muts in mutations_by_type.items():
                if mutation_type in mutation_symbols:
                    symbol_info = mutation_symbols[mutation_type]
                    positions = [mut['position'] for mut in muts]
                    # Enhanced hover text with UMI ID and sample info
                    hover_text = []
                    for mut in muts:
                        hover_text.append(f"Sample: {sample_name}<br>UMI: {consensus_id}<br>Position: {mut['position']}<br>Mutation: {mut['ref_base']}→{mut['cons_base']}<br>Type: {mut['type']}")
                    
                    fig.add_trace(go.Scatter(
                        x=positions,
                        y=[y_offset] * len(positions),
                        mode='markers',
                        marker=dict(
                            color=symbol_info['color'],
                            size=12,
                            symbol=symbol_info['symbol'],
                            line=dict(width=2, color='black')
                        ),
                        name=mutation_type,  # Only show mutation type in legend
                        showlegend=False,  # Hide from legend to avoid duplicates
                        hovertemplate='<br>'.join(hover_text) + '<extra></extra>'
                    ))
                else:
                    # Default symbol for unknown mutation types
                    positions = [mut['position'] for mut in muts]
                    hover_text = []
                    for mut in muts:
                        hover_text.append(f"Sample: {sample_name}<br>UMI: {consensus_id}<br>Position: {mut['position']}<br>Mutation: {mut['ref_base']}→{mut['cons_base']}<br>Type: {mut['type']}")
                    
                    fig.add_trace(go.Scatter(
                        x=positions,
                        y=[y_offset] * len(positions),
                        mode='markers',
                        marker=dict(
                            color='red',
                            size=12,
                            symbol='x',
                            line=dict(width=2, color='black')
                        ),
                        name=mutation_type,  # Only show mutation type in legend
                        showlegend=False,  # Hide from legend to avoid duplicates
                        hovertemplate='<br>'.join(hover_text) + '<extra></extra>'
                    ))
        
        y_offset += 0.15  # Much tighter spacing between reads
    
    # Add position markers every 500 bases for reference
    position_markers = list(range(0, ref_length, 500))
    if position_markers:
        fig.add_trace(go.Scatter(
            x=position_markers,
            y=[-0.5] * len(position_markers),
            mode='markers',
            marker=dict(color='gray', size=4, symbol='line-ns'),
            name='Position Markers',
            showlegend=False,
            hovertemplate='Position: %{x}<extra></extra>'
        ))
    
    fig.update_layout(
        title="Consensus Reads Alignment with Mutations",
        xaxis_title="Genome Position (bp)",
        yaxis_title="Reads",
        height=max(500, len(alignment_df) * 60 + 300),
        showlegend=True,
        hovermode='closest',
        xaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'
        ),
        yaxis=dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray'
        )
    )
    
    return html.Div([
        html.H4("🧬 Alignment Visualizer"),
        html.P("Interactive alignment showing consensus reads aligned to reference sequence with mutations highlighted in red."),
        dcc.Graph(figure=fig),
        html.H5("Alignment Details (Sortable)"),
        dash_table.DataTable(
            data=alignment_df.to_dict('records'),
            columns=[
                {"name": "Sample", "id": "sample", "type": "text"},
                {"name": "Consensus ID", "id": "consensus_id", "type": "text"},
                {"name": "Mutations Count", "id": "mutations_count", "type": "numeric"},
                {"name": "Alignment Score", "id": "alignment_score", "type": "numeric", "format": {"specifier": ".1f"}},
                {"name": "Matches", "id": "matches", "type": "numeric"},
                {"name": "Mismatches", "id": "mismatches", "type": "numeric"},
                {"name": "Gaps", "id": "gaps", "type": "numeric"},
                {"name": "Mutation Details", "id": "mutation_details", "type": "text"}
            ],
            sort_action="native",
            sort_mode="multi",
            filter_action="native",
            page_action="native",
            page_current=0,
            page_size=20,
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ]
        )
    ])


def create_raw_data_tab():
    """Create the raw data tab"""
    global analysis_data
    
    if not analysis_data:
        return html.Div("Analysis data not available.")
    
    # Display raw JSON data
    raw_data = json.dumps(analysis_data, indent=2, default=str)
    
    return html.Div([
        html.H4("Raw Analysis Data"),
        html.Pre(raw_data, style={'whiteSpace': 'pre-wrap', 'wordBreak': 'break-all'})
    ])




def main():
    """Main function to run the web application"""
    import argparse
    
    parser = argparse.ArgumentParser(description="BELT Mutation Analysis Visualizer")
    parser.add_argument("--port", type=int, default=8050, help="Port to run the web app on")
    parser.add_argument("--host", default="127.0.0.1", help="Host to run the web app on")
    parser.add_argument("--debug", action="store_true", help="Run in debug mode")
    
    args = parser.parse_args()
    
    logger.info(f"Starting BELT Mutation Analysis Visualizer on {args.host}:{args.port}")
    app.run(debug=args.debug, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
