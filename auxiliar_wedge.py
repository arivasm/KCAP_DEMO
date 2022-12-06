import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from pyDatalog import pyDatalog
from pyDatalog.pyDatalog import assert_fact, load, ask

from SPARQLWrapper import SPARQLWrapper, JSON
from rdflib import Graph
import json

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.patches import FancyArrowPatch

from ipycytoscape import *
import ipywidgets as widgets
import bqplot
from bqplot import pyplot as bqplt
from bqplot import Tooltip
from scipy.stats import wilcoxon, kruskal


def build_query_clarify(input_cui, type_ddi):
    input_cui_uri = ','.join(['<http://clarify2020.eu/entity/' + cui + '>' for cui in input_cui])
    query = """
    select distinct ?EffectorDrugLabel ?AffectedDrugLabel ?Effect ?Impact ?precipitantDrug ?objectDrug
        where {""" + type_ddi + """        
        ?s <http://clarify2020.eu/vocab/effect_cui> ?o . 
        ?o <http://clarify2020.eu/vocab/annLabel> ?Effect . 
        ?s <http://clarify2020.eu/vocab/impact> ?Impact .
        ?s <http://clarify2020.eu/vocab/precipitant_drug_cui> ?precipitantDrug .
        ?s <http://clarify2020.eu/vocab/object_drug_cui> ?objectDrug .
        ?precipitantDrug <http://clarify2020.eu/vocab/annLabel> ?EffectorDrugLabel.
        ?objectDrug <http://clarify2020.eu/vocab/annLabel> ?AffectedDrugLabel.
    FILTER (?precipitantDrug in (""" + input_cui_uri + """ ) && ?objectDrug in (""" + input_cui_uri + """))
    }"""
    return query


def store_pharmacokinetic_ddi(effect):
    if effect in ['Excretion_rate', 'Excretory_function', 'Excretion']:
        effect = 'excretion'
    elif effect in ['Process_of_absorption', 'Absorption']:
        effect = 'absorption'
    elif effect in ['Serum_concentration', 'Serum_concentration_of', 'Serum_level', 'Serum_globulin_level', 'Metabolite', 'Active_metabolites']:
        effect = 'serum_concentration'
    elif effect in ['Metabolism']:
        effect = 'metabolism'
    # else:
    #    return 'pharmacodynamic'
    return effect


def rename_impact(impact):
    if impact in ['Increase', 'Higher', 'Worsening']:
        return 'increase'
    return 'decrease'


def query_result_clarify(query, endpoint):
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'EffectorDrugLabel': [], 'AffectedDrugLabel': [], 'Effect': [], 'Impact': [], 'precipitantDrug': [],
          'objectDrug': []}
    for r in results['results']['bindings']:
        effect = r['Effect']['value']
        effect = store_pharmacokinetic_ddi(effect)
        if effect == 'pharmacodynamic':
            continue
        dd['Effect'].append(effect.lower())
        impact = r['Impact']['value'].replace('http://clarify2020.eu/entity/', '')
        impact = rename_impact(impact)
        dd['Impact'].append(impact)
        dd['EffectorDrugLabel'].append(r['EffectorDrugLabel']['value'].lower())
        dd['AffectedDrugLabel'].append(r['AffectedDrugLabel']['value'].lower())
        dd['precipitantDrug'].append(r['precipitantDrug']['value'].replace('http://clarify2020.eu/entity/', ''))
        dd['objectDrug'].append(r['objectDrug']['value'].replace('http://clarify2020.eu/entity/', ''))

    set_DDIs = pd.DataFrame(dd)

    return set_DDIs


def query_result_symmetric_clarify(query, endpoint):
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    dd = {'EffectorDrugLabel': [], 'AffectedDrugLabel': [], 'Effect': [], 'Impact': [], 'precipitantDrug': [],
          'objectDrug': []}
    for r in results['results']['bindings']:
        effect = r['Effect']['value']
        effect = store_pharmacokinetic_ddi(effect)
        if effect == 'pharmacodynamic':
            continue
        dd['Effect'].append(effect.lower())
        impact = r['Impact']['value'].replace('http://clarify2020.eu/entity/', '')
        impact = rename_impact(impact)
        dd['Impact'].append(impact)
        dd['EffectorDrugLabel'].append(r['AffectedDrugLabel']['value'].lower())
        dd['AffectedDrugLabel'].append(r['EffectorDrugLabel']['value'].lower())
        dd['precipitantDrug'].append(r['objectDrug']['value'].replace('http://clarify2020.eu/entity/', ''))
        dd['objectDrug'].append(r['precipitantDrug']['value'].replace('http://clarify2020.eu/entity/', ''))

    set_DDIs = pd.DataFrame(dd)
    #lowerify_cols = [col for col in set_DDIs if col not in ['precipitantDrug', 'objectDrug']]
    #set_DDIs[lowerify_cols] = set_DDIs[lowerify_cols].apply(lambda x: x.astype(str).str.lower(), axis=1)
    return set_DDIs


def combine_col(corpus, cols):
    # corpus = corpus.apply(lambda x: x.astype(str).str.lower())
    name = '_'.join(cols)
    corpus[name] = corpus[cols].apply(lambda x: '_'.join(x.values.astype(str)), axis=1)
    return corpus


def get_drug_label_by_category(drugs_cui, set_DDIs):
    d_label = set(set_DDIs.loc[set_DDIs.precipitantDrug.isin(drugs_cui)].EffectorDrugLabel.unique())
    d_label.update(set_DDIs.loc[set_DDIs.objectDrug.isin(drugs_cui)].AffectedDrugLabel.unique())
    # drugs_cui = [c.lower() for c in drugs_cui]
    # print(drugs_cui)
    # d_label = set(set_DDIs.loc[set_DDIs.EffectorDrugLabel.isin(drugs_cui)].EffectorDrugLabel.unique())
    # d_label.update(set_DDIs.loc[set_DDIs.AffectedDrugLabel.isin(drugs_cui)].AffectedDrugLabel.unique())
    return d_label


def extract_ddi(onco_drugs, non_onco_drugs, endpoint):
    input_cui = onco_drugs + non_onco_drugs

    # ===== Endpoint 'clarify2020' =====
    asymmetric_ddi = """    ?s a <http://clarify2020.eu/vocab/DrugDrugInteraction> .
                    """
    symmetric_ddi = """     ?sim a <http://clarify2020.eu/vocab/SymmetricDrugDrugInteraction> .
                           ?sim <http://www.w3.org/2000/01/rdf-schema#subClassOf> ?s.
                   """
    # ===== Endpoint 'covid-19' =====
    # asymmetric_ddi = """    ?s a <http://covid-19.tib.eu/vocab/DrugDrugInteraction> .
    #                  """
    # symmetric_ddi = """     ?sim a <http://covid-19.tib.eu/vocab/SymmetricDrugDrugInteraction> .
    #                         ?sim <http://www.w3.org/2000/01/rdf-schema#subClassOf> ?s.
    #                 """
    # prefix = 'http://covid-19.tib.eu/Annotation/'

    query = build_query_clarify(input_cui, asymmetric_ddi)
    print(query)
    set_DDIs = query_result_clarify(query, endpoint)
    query = build_query_clarify(input_cui, symmetric_ddi)
    corpus_symmetric = query_result_symmetric_clarify(query, endpoint)

    set_DDIs = combine_col(set_DDIs, ['Effect', 'Impact'])
    corpus_symmetric = combine_col(corpus_symmetric, ['Effect', 'Impact'])
    adverse_event = corpus_symmetric.Effect_Impact.unique()
    union = pd.concat([set_DDIs, corpus_symmetric])

    set_dsd_label = get_drug_label_by_category(onco_drugs, union)
    comorbidity_drug = get_drug_label_by_category(non_onco_drugs, union)
    set_DDIs = set_DDIs[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    return adverse_event, union, set_dsd_label, comorbidity_drug, set_DDIs


def load_data(file):
    #onco_drugs = file["Input"]["OncologicalDrugs"]
    #non_onco_drugs = file["Input"]["Non_OncologicalDrugs"]
    onco_drugs = file["oncological_drug"]
    non_onco_drugs = file["non_oncological_drug"]
    return extract_ddi(onco_drugs, non_onco_drugs, "https://labs.tib.eu/sdm/clarify-kg-8-1/sparql")


pyDatalog.create_terms('rdf_star_triple, inferred_rdf_star_triple, wedge, A, B, C, T, T2')

"""
def build_datalog_model(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T) & (T._in(ddiTypeToxicity))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeToxicity)) & (T2._in(ddiTypeToxicity)) & (A != C)
    wedge(A, B, C, T, T2) <= inferred_rdf_star_triple(A, B, T) & inferred_rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeToxicity)) & (T2._in(ddiTypeToxicity)) & (A != C)

    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T) & (T._in(ddiTypeEffectiveness))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeEffectiveness)) & (T2._in(ddiTypeEffectiveness)) & (A != C)
    wedge(A, B, C, T, T2) <= inferred_rdf_star_triple(A, B, T) & inferred_rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeEffectiveness)) & (T2._in(ddiTypeEffectiveness)) & (A != C)


def compute_wedge_datalog(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    wedge(A, B, C, T, T2) <= rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (T2._in(ddiTypeToxicity)) & (
        T._in(ddiTypeToxicity)) & (A != C)
    wedge(A, B, C, T, T2) <= rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (T2._in(ddiTypeEffectiveness)) & (
        T._in(ddiTypeEffectiveness)) & (A != C)
"""


def build_datalog_model(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T) # & (T._in(ddiTypeToxicity))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeToxicity)) & (T2._in(ddiTypeToxicity)) & (A != C)

    inferred_rdf_star_triple(A, B, T) <= rdf_star_triple(A, B, T) # & (T._in(ddiTypeEffectiveness))
    inferred_rdf_star_triple(A, C, T2) <= inferred_rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (
        T._in(ddiTypeEffectiveness)) & (T2._in(ddiTypeEffectiveness)) & (A != C)

    wedge(A, B, C, T, T2) <= inferred_rdf_star_triple(A, B, T) & inferred_rdf_star_triple(B, C, T2) & (A != C)


def compute_wedge_datalog(union):
    pyDatalog.clear()
    for d in union.values:
        # Extensional Database
        assert_fact('rdf_star_triple', d[0], d[1], d[2])
    # Intentional Database
    wedge(A, B, C, T, T2) <= rdf_star_triple(A, B, T) & rdf_star_triple(B, C, T2) & (A != C)

    
def get_indirect_ddi(indirect_ddi, dsd, write=False):
    derived_ddi = []
    deduced_ddi = inferred_rdf_star_triple(C, dsd, T)
    for i in range(len(deduced_ddi)):
        x = {'EffectorDrugLabel': [deduced_ddi[i][0]], 'AffectedDrugLabel': dsd,
             'Effect_Impact': deduced_ddi[i][1]}  # + '_derived'
        indirect_ddi = pd.concat([indirect_ddi, pd.DataFrame(data=x)])
        if write:
            impact = deduced_ddi[i][1].split('_')[-1]
            effect = deduced_ddi[i][1].split('_')[:-1]
            effect = '_'.join([l for l in effect])
            # print(effect, impact)
            derived_ddi.append(deduced_ddi[i][0] + ' ' + impact + ' ' + effect + ' of ' + dsd + ' (derived)')

    return indirect_ddi, derived_ddi


def get_indirect_ddi_treatment(set_dsd_label, write):
    indirect_ddi = pd.DataFrame(columns=['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact'])
    text_derived_ddi = []
    for dsd in set_dsd_label:
        indirect_ddi, derived_ddi = get_indirect_ddi(indirect_ddi, dsd, write)
        text_derived_ddi = text_derived_ddi + derived_ddi
    return indirect_ddi, text_derived_ddi


# # Add color to edges
def add_color(set_DDIs):
    ColorLegend = {}

    col = colors.cnames
    color = list(col.values())

    # --------------remove colors-------------
    index = color.index('#FFEBCD')
    color.pop(index)
    index = color.index('#7FFF00')
    color.pop(index)
    index = color.index('#FFF8DC')
    color.pop(index)
    index = color.index('#A9A9A9')
    color.pop(index)
    # --------------remove colors-------------

    effect_impact = list(pd.value_counts(set_DDIs.Effect_Impact).index)
    # print(effect_impact)
    for i in range(len(effect_impact)):
        set_DDIs.loc[set_DDIs.Effect_Impact == effect_impact[i], "edge_color"] = color[i + 8]
        set_DDIs.loc[set_DDIs.Effect_Impact == effect_impact[i] + '_derived', "edge_color"] = color[i + 8]
        set_DDIs.loc[set_DDIs.Effect_Impact == effect_impact[i].replace('_derived', ''), "edge_color"] = color[i + 8]

        ColorLegend[effect_impact[i]] = color[i + 8]
        # print(effect_impact[i], color[i + 8])

    set_DDIs.reset_index(inplace=True)
    set_DDIs = set_DDIs.drop(columns=['index'])
    return ColorLegend, set_DDIs


# # Visualize graph
def add_di_edge_to_graph(set_ddi, g):
    # Add edges and edge attributes
    for i, elrow in set_ddi.iterrows():
        g.add_edge(elrow[0], elrow[1], ddi_type=elrow[2], ddi_color=elrow[3])  # , create_using=nx.MultiDiGraph()
    # g.add_node('enalapril')
    return g


def get_multiple_edge(g):
    multiple_edge = []
    for e in g.edges:
        multiple_edge.append(0.12 * (e[2]))
    return multiple_edge


def get_node_size(g):
    list_node = list(g.nodes)
    d2_size = []
    for n2 in list_node:
        d2_size.append(150 + (g.out_degree(n2) * 50))
    return d2_size


# # Color nodes
def get_node_color(g, label_dsd, color_mark, color_basic):
    list_node = list(g.nodes)
    # print('check:', list_node)
    n_color = [color_basic] * g.number_of_nodes()
    for d in label_dsd:
        if d in list_node:
            index = list_node.index(d)
            n_color[index] = color_mark
    return n_color


def get_effect_impact(edge):
    l = edge.split('_')
    impact = l[len(l) - 1]
    effect = '_'.join(l[:len(l) - 1])
    return effect, impact


def plot_graph(g, n_color, d2_size, multiple_edge, ColorLegend, plot_name, dim, mechanism, adverse_event,
               graph_enriched=True):
    increased_ddi = []
    fig = plt.figure(figsize=dim)  # 7,4
    ax = plt.subplot(1, 1, 1)
    pos = nx.circular_layout(g)
    nx.draw_networkx_nodes(g, pos=pos, node_color=n_color, node_size=d2_size)  # 'skyblue'  , label=True

    list_annotation = []
    list_ddi = []
    i = 0
    for edge in g.edges(data=True):
        dash = False
        ls_dash = None
        label_ddi = edge[2]['ddi_type']

        for pair in g.in_edges(edge[0]):
            if dash: break
            for e in g.edges(data=True):
                if e[0] == pair[0] and e[1] == pair[1] and e[2]['ddi_type'] in mechanism:
                    dash = True
                    ls_dash = 'dotted'  # dashed, dashdot, dotted
                    label_ddi = label_ddi + '**'
                    if not edge[2]['ddi_type'].endswith('_derived'):
                        effect, impact = get_effect_impact(e[2]['ddi_type'])
                        antecedent = ' is increased because ' + e[0] + ' ' + impact + ' ' + effect + ' of ' + e[1]
                        effect, impact = get_effect_impact(edge[2]['ddi_type'])
                        increased_ddi.append(
                            'The interaction: ' + edge[0] + ' ' + impact + ' ' + effect + ' of ' + edge[1] + antecedent)

                    break
        if graph_enriched:
            if edge[2]['ddi_type'] in mechanism:
                a1 = ax.annotate("", xy=pos[edge[1]], xytext=pos[edge[0]],
                                 arrowprops=dict(arrowstyle="fancy,head_width=0.2,tail_width=0.6",
                                                 color=edge[2]['ddi_color'],
                                                 ls=ls_dash, shrinkA=10, shrinkB=10, fc="1.0", lw=2.0,
                                                 connectionstyle="arc3,rad=" + str(multiple_edge[i]), ),
                                 label=label_ddi)  # fc="0.0",  +str('+')

            elif edge[2]['ddi_type'] in adverse_event:
                a1 = ax.annotate("", xy=pos[edge[1]], xytext=pos[edge[0]],
                                 arrowprops=dict(arrowstyle='-', color=edge[2]['ddi_color'],
                                                 shrinkA=5, shrinkB=5, lw=0.5,
                                                 connectionstyle="arc3,rad=" + str(multiple_edge[i]), ),
                                 label=label_ddi)

            elif not edge[2]['ddi_type'].endswith('_derived'):
                a1 = ax.annotate("", xy=pos[edge[1]], xytext=pos[edge[0]],
                                 arrowprops=dict(arrowstyle='-|>, head_length=0.6, head_width=0.3',
                                                 color=edge[2]['ddi_color'],
                                                 ls=ls_dash, shrinkA=3, shrinkB=3, lw=2.0,
                                                 connectionstyle="arc3,rad=" + str(multiple_edge[i]), ),
                                 label=label_ddi)

            else:  # edge[2]['ddi_type'].endswith('_derived'):
                a1 = ax.annotate("", xy=pos[edge[1]], xytext=pos[edge[0]],
                                 arrowprops=dict(arrowstyle='-|>, head_length=0.6, head_width=0.3',
                                                 color=edge[2]['ddi_color'],
                                                 shrinkA=5, shrinkB=5, lw=0.5,
                                                 connectionstyle="arc3,rad=" + str(multiple_edge[i]), ),
                                 label=label_ddi)
        else:
            a1 = ax.annotate("", xy=pos[edge[1]], xytext=pos[edge[0]],
                             arrowprops=dict(arrowstyle='-|>, head_length=0.6, head_width=0.3',
                                             color=edge[2]['ddi_color'],
                                             shrinkA=5, shrinkB=5, lw=0.5,
                                             connectionstyle="arc3,rad=" + str(multiple_edge[i]), ),
                             label=edge[2]['ddi_type'])

        i += 1
        if edge[2]['ddi_type'] in list_ddi:
            continue
        list_annotation.append(a1)
        list_ddi.append(edge[2]['ddi_type'])

    # sort the legend by length of the string
    sorted_list = sorted(list_ddi, key=len)
    copy_list_annotation = list_annotation.copy()
    for i in range(len(list_ddi)):
        j = sorted_list.index(list_ddi[i])
        list_annotation[j] = copy_list_annotation[i]

    nx.draw_networkx_labels(g, pos=pos, font_size=14, alpha=1.0)  # ,font_color='r'  , font_weight="bold"

    # ax.margins(0.1, 0.1)
    # ax.set_xmargin(0.35)
    plt.axis('off')
    plt.title('Drug-Drug Interactions', fontsize=12, ha='center')

    ax.legend(handles=list(list_annotation), handler_map={type(a1): AnnotationHandler(5)},
              ncol=2, bbox_to_anchor=(0.5, 0), loc='upper center').get_frame().set_alpha(0.0)
    ax.set_xmargin(0.35)
    plt.tight_layout()
    plt.savefig('output/' + plot_name + '.png', format='png', dpi=300)  # , format='png', dpi = 300   .pdf
    # plt.show()
    return increased_ddi, fig


def visualise_treatment(diff_df, sd, plot_name, dim, mechanism, adverse_event, plot_treatment=True,
                        graph_enriched=True):
    g_ddi = diff_df.copy()
    increased_ddi = []
    ColorLegend, union = add_color(g_ddi)  # set_unsymmetric_ddi
    g = nx.MultiDiGraph()  # DiGraph

    g = add_di_edge_to_graph(union, g)
    if plot_treatment:
        multiple_edge = get_multiple_edge(g)
        if sd == 'all_drug':
            n_color = 'skyblue'
        else:
            n_color = get_node_color(g, sd, 'red', 'skyblue')
        d2_size = get_node_size(g)
        increased_ddi = plot_graph(g, n_color, d2_size, multiple_edge, ColorLegend, plot_name, dim, mechanism,
                                   adverse_event, graph_enriched)
    return g, increased_ddi


class AnnotationHandler(HandlerLine2D):
    def __init__(self, ms, *args, **kwargs):
        self.ms = ms
        HandlerLine2D.__init__(self, *args, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        orig_handle.arrowprops['shrinkA'] = 5
        orig_handle.arrowprops['shrinkB'] = 5
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent, width, height, fontsize)
        ydata = ((height - ydescent) / 2.) * np.ones(len(xdata), float)
        legline = FancyArrowPatch(posA=(xdata[0], ydata[0]), posB=(xdata[-1], ydata[-1]),
                                  mutation_scale=self.ms, **orig_handle.arrowprops)
        legline.set_transform(trans)
        return legline,


def create_json_to_cytoscape(union, k):
    graph_json = dict()
    graph_json['nodes'] = []
    graph_json['edges'] = []
    drug_id = dict()
    id_x = 0
    for i in range(union.shape[0]):
        precipitant = union.iloc[i]['EffectorDrugLabel']
        object_d = union.iloc[i]['AffectedDrugLabel']
        ddi = union.iloc[i]['Effect_Impact']
        edge = dict()
        edge['data'] = dict()

        if precipitant in drug_id.keys():
            edge['data']['id'] = id_x
            edge['data']['source'] = drug_id[precipitant]
            edge['data']['Effect_Impact'] = ddi
            id_x += 1
        else:
            node = dict()
            node['data'] = dict()
            drug_id[precipitant] = id_x
            node['data']['id'] = id_x
            node['data']['name'] = precipitant
            edge['data']['id'] = id_x + 1
            edge['data']['source'] = id_x
            edge['data']['Effect_Impact'] = ddi
            graph_json['nodes'].append(node)
            id_x += 2
        if object_d in drug_id.keys():
            edge['data']['target'] = drug_id[object_d]
        else:
            node = dict()
            node['data'] = dict()
            drug_id[object_d] = id_x
            node['data']['id'] = id_x
            node['data']['name'] = object_d
            edge['data']['target'] = id_x
            graph_json['nodes'].append(node)
            id_x += 1
            if object_d == k:
                node['classes'] = 'red'  # Single class

        graph_json['edges'].append(edge)

    return graph_json


# # Whole Graph enriched
def get_graph_enriched(plot_ddi, comorbidity_drug, set_dsd_label, adverse_event, set_DDIs):
    build_datalog_model(plot_ddi)
    indirect_ddi, text_derived_ddi = get_indirect_ddi_treatment(comorbidity_drug.union(set_dsd_label), write=True)
    mechanism = ddiTypeEffectiveness + ddiTypeToxicity
    g_i, increased_ddi = visualise_treatment(plot_ddi, list(set_dsd_label), 'Graph_initial', (7, 5), mechanism, adverse_event,
                        plot_treatment=True, graph_enriched=False)
    graph_ddi = pd.concat([set_DDIs, indirect_ddi])
    graph_ddi.drop_duplicates(keep='first', inplace=True)
    g_j, increased_ddi = visualise_treatment(graph_ddi, list(set_dsd_label), 'Graph_enriched', (7, 5), mechanism,
                                           adverse_event)
    
    # centrality
    #deg_centrality = nx.degree_centrality(g)
    #centrality = np.fromiter(deg_centrality.values(), float)
    #print(deg_centrality, 'centrality:', centrality)
    #print(nx.clustering(g_j))
    #G2 = nx.Graph(g_j)
    #print(nx.average_clustering(G2))
    #print('density: ', nx.density(g))
    #print('closeness_centrality: ',nx.closeness_centrality(g))
    #print('degree: ', nx.degree(g))
    #print('info: ', nx.info(g))
    #for i in nx.strongly_connected_components(g):
    #    print('strongly_connected_components: ', i)
    #print('connected_components: ', nx.number_strongly_connected_components(g))
    #print('average_neighbor_degree: ', nx.average_neighbor_degree(g))
    #print('simple_cycles: ', list(nx.simple_cycles(g)))
    #print('diameter: ', nx.diameter(g.to_undirected())) # it's not a good metric for this problem!
    comparision = pd.DataFrame(index=['Density', 'Num. edges', 'Num. cycles', 'Avg. clustering coefficient', 'Connected Component'],
                           columns=['Graph Initial', 'Graph Enriched'])
    comparision.at['Density','Graph Initial']=nx.density(g_i)
    comparision.at['Num. edges','Graph Initial']=nx.number_of_edges(g_i)
    G2 = nx.Graph(g_i)
    comparision.at['Avg. clustering coefficient','Graph Initial']=nx.average_clustering(G2)
    comparision.at['Num. cycles','Graph Initial']=len(list(nx.simple_cycles(g_i)))
    comparision.at['Connected Component','Graph Initial']=nx.number_strongly_connected_components(g_i)
    
    comparision.at['Density','Graph Enriched']=nx.density(g_j)    
    comparision.at['Num. edges','Graph Enriched']=nx.number_of_edges(g_j)
    G2 = nx.Graph(g_j)
    comparision.at['Avg. clustering coefficient','Graph Enriched']=nx.average_clustering(G2)   
    comparision.at['Num. cycles','Graph Enriched']=len(list(nx.simple_cycles(g_j)))
    comparision.at['Connected Component','Graph Enriched']=nx.number_strongly_connected_components(g_j) 
    #display(comparision)
    return comparision, graph_ddi, text_derived_ddi


def get_all_wedge(union):
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    build_datalog_model(plot_ddi)
    w = wedge(A, B, C, T, T2)
    return w


def computing_wedge(set_drug_label):
    dict_wedge = dict()
    dict_frequency = dict()
    for d in set_drug_label:
        w = wedge(A, d, C, T, T2)
        indirect_ddi = pd.DataFrame(columns=['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact'])
        for i in range(len(w)):
            x = {'EffectorDrugLabel': [w[i][0], d], 'AffectedDrugLabel': [d, w[i][1]],
                 'Effect_Impact': [w[i][2], w[i][3]]}
            indirect_ddi = pd.concat([indirect_ddi, pd.DataFrame(data=x)])

        indirect_ddi.drop_duplicates(keep='first', inplace=True)
        dict_wedge[d] = indirect_ddi
        dict_frequency[d] = len(w)
    return dict_wedge, dict_frequency


ddiTypeToxicity = ["serum_concentration_increase", "metabolism_decrease", "absorption_increase", "excretion_decrease"]
ddiTypeEffectiveness = ["serum_concentration_decrease", "metabolism_increase", "absorption_decrease",
                        "excretion_increase"]


# adverse_event, union, set_dsd_label, comorbidity_drug, set_DDIs = load_data(file)


def capture_knowledge(adverse_event, union, set_dsd_label, comorbidity_drug, set_DDIs):
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    network_analysis, graph_ddi, text_derived_ddi = get_graph_enriched(plot_ddi, comorbidity_drug, set_dsd_label,
                                                                    adverse_event, set_DDIs)  #plot_ddi
    #print('Num.DDIs Initial: ', plot_ddi.shape[0], 'Num.DDIs Enriched: ', graph_ddi.shape[0])
    return network_analysis, graph_ddi

def discovering_knowledge(adverse_event, union, set_dsd_label, comorbidity_drug):
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    build_datalog_model(plot_ddi)
    dict_wedge, dict_frequency = computing_wedge(comorbidity_drug.union(set_dsd_label))
    dict_frequency = dict(sorted(dict_frequency.items(), key=lambda item: item[1], reverse=True))

    dict_graph_json = dict()
    for k, v in dict_frequency.items():
        if v > 0:
            #visualise_treatment(dict_wedge[k], {k}, 'middle_vertex' + k, (7, 5), ddiTypeToxicity, adverse_event,
            #                    plot_treatment=True, graph_enriched=False)
            graph_json = create_json_to_cytoscape(dict_wedge[k], k)
            dict_graph_json[k] = graph_json
    return dict_graph_json, dict_frequency


def evaluation_without_deduction(union, set_dsd_label, comorbidity_drug):
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    compute_wedge_datalog(plot_ddi)
    dict_wedge, dict_frequency = computing_wedge(comorbidity_drug.union(set_dsd_label))
    dict_frequency = dict(sorted(dict_frequency.items(), key=lambda item: item[1], reverse=True))
    return dict_frequency


def create_graph_cytoscape(middle_vertex):
    # load a style dictionary
    with open("styles.json") as fi:
        s = json.load(fi)
    # Create the cytoscape graph widget
    cytoscapeobj = CytoscapeWidget()
    cytoscapeobj.graph.add_graph_from_json(middle_vertex, directed=True, multiple_edges=True)  # , directed=True, input_data['elements']
    
    cytoscapeobj.set_style(s)
    cytoscapeobj.set_layout(name='breadthfirst', animate=True, nodeSpacing = 5)  # concentric,  breadthfirst
    return cytoscapeobj


def distribution_wedge(dict_frequency):
    plt.bar(*zip(*dict_frequency.items()))
    plt.title('Wedge frequency distribution')
    plt.ylabel('Wedge absolute frequency')
    plt.xlabel('Middle-Vertex')
    plt.show()


def comparision_distribution_wedge(df1):
    fig = bqplt.figure(title="Wedge frequency distribution", fig_margin={'top':40, 'bottom':40, 'left':60, 'right':0},
                     legend_location="top-right")

    bar_chart  = bqplt.bar(x = list(df1.index), y= [df1["Eval1"], df1["Eval2"]],
                         labels = ["Capturing Knowledge", "Baseline"], display_legend=True)

    bar_chart.type = "grouped"
    bar_chart.colors = ["tomato", "DeepSkyBlue"]
    bar_chart.tooltip = Tooltip(fields=["x", "y"], labels=["Middle-Vertex", "Wedge absolute frequency"])
    bqplt.xlabel("Middle-Vertex")
    bqplt.ylabel("Wedge absolute frequency")
    bqplt.show()


# Wilcoxon signed-rank test
def wilcoxon_test(experiment_result):
    data1 = experiment_result.Eval1.values
    data2 = experiment_result.Eval2.values
    if set(data1)==set(data2):
        print('Same distributions (fail to reject H0)', '\nStatistics=0.0000, p=1.0000000000')
        return
    # compare samples
    stat, p = wilcoxon(data1, data2) #, mode='approx'
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Same distributions (fail to reject H0)')
    else:
        print('Different distributions (reject H0)')

    print('Statistics=%.4f, p=%.10f' % (stat, p))


def kruskal_test(experiment_result):
    data1 = experiment_result.Eval1.values
    data2 = experiment_result.Eval2.values
    # compare samples
    stat, p = kruskal(data1, data2)
    # interpret
    alpha = 0.05
    if p > alpha:
        print('Same distributions (fail to reject H0)')
    else:
        print('Different distributions (reject H0)')

    print('Statistics=%.4f, p=%.10f' % (stat, p))
    

# # Graph with minimal set PCD withdrawn
def get_minimal_set(comorbidity_drug, set_dsd_label, union, adverse_event):
    plot_ddi = union[['EffectorDrugLabel', 'AffectedDrugLabel', 'Effect_Impact']]
    plot_ddi.drop_duplicates(keep='first', inplace=True)
    minimal_set_d = {}
    mechanism = ddiTypeEffectiveness + ddiTypeToxicity
    for comorb_d in comorbidity_drug:

        minimal_set = plot_ddi.loc[plot_ddi.EffectorDrugLabel != comorb_d]
        minimal_set = minimal_set.loc[minimal_set.AffectedDrugLabel != comorb_d]
        if minimal_set.empty:
            continue
        build_datalog_model(minimal_set)
        indirect_ddi, derived_ddi = get_indirect_ddi_treatment(comorbidity_drug.union(set_dsd_label), False)
        graph_ddi = pd.concat([minimal_set, indirect_ddi])
        graph_ddi.drop_duplicates(keep='first', inplace=True)
        visualise_treatment(graph_ddi, list(set_dsd_label), 'name', (6, 5), mechanism, adverse_event,
                            plot_treatment=False)
        # print(comorb_d, graph_ddi.shape[0])
        if not minimal_set_d:
            minimal_set_d[comorb_d] = {}
            minimal_set_d[comorb_d]['n_ddi'] = graph_ddi.shape[0]
            minimal_set_d[comorb_d]['graph'] = graph_ddi
        elif list(minimal_set_d.values())[0]['n_ddi'] > graph_ddi.shape[0]:
            minimal_set_d = {}
            minimal_set_d[comorb_d] = {}
            minimal_set_d[comorb_d]['n_ddi'] = graph_ddi.shape[0]
            minimal_set_d[comorb_d]['graph'] = graph_ddi
        elif list(minimal_set_d.values())[0]['n_ddi'] == graph_ddi.shape[0]:
            minimal_set_d[comorb_d] = {}
            minimal_set_d[comorb_d]['n_ddi'] = graph_ddi.shape[0]
            minimal_set_d[comorb_d]['graph'] = graph_ddi

    if len(minimal_set_d) > 1:
        max_ddi_mechanism = pd.DataFrame()
        result = {}
        for key, value in minimal_set_d.items():
            ddi_mechanism = plot_ddi.loc[((plot_ddi.EffectorDrugLabel == key) | (plot_ddi.AffectedDrugLabel == key)) &
                                         (plot_ddi.Effect_Impact.isin(mechanism))]
            print(key, ddi_mechanism)
            if ddi_mechanism.shape[0] > max_ddi_mechanism.shape[0]:
                max_ddi_mechanism = ddi_mechanism
                result = {key: value}
                drug_remove = key
                mechanism_flag = True
        if mechanism_flag:
            graph_ddi = result[drug_remove]['graph']
        else:
            graph_ddi = minimal_set_d.keys()[0]['graph']
            drug_remove = minimal_set_d.keys()[0]
        print('minimal set of PCD: ', drug_remove, 'Num.DDIs: ', graph_ddi.shape[0])
        visualise_treatment(graph_ddi, list(set_dsd_label), 'minimal_set', (6, 5), mechanism, adverse_event)

    elif len(minimal_set_d) == 0:
        return 0, "Removing any comorbidity drugs there is no interaction in treatment."
    else:
        drug_remove = list(minimal_set_d.keys())[0]
        graph_ddi = minimal_set_d[drug_remove]['graph']
        print('minimal set of PCD: ', drug_remove, 'Num.DDIs: ', graph_ddi.shape[0])
        visualise_treatment(graph_ddi, list(set_dsd_label), 'minimal_set', (6, 5), mechanism, adverse_event)
    return graph_ddi, drug_remove

