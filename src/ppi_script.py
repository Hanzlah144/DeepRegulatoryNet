import os
import pandas as pd
import requests
import networkx as nx
from io import StringIO
from pyvis.network import Network
from pathlib import Path
import logging
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
ppi_logger = logging.getLogger(__name__)

try:
    from IPython.display import display, HTML
    boolean_track = True
except ImportError:
    boolean_track = False
    ppi_logger.warning("[WARN]  IPython.display not available")


def gene_path(gene_path):
    logger = logging.getLogger()
    try:
        gene_df = pd.read_csv(gene_path)
        if 'Gene' not in gene_df.columns:
            raise ValueError("CSV file should contain a 'Gene' column")
        genes = gene_df['Gene'].dropna().tolist()
        return [gene for gene in genes if isinstance(gene, str)]
    except Exception as e:
        logger.error("[ERROR]  Error reading gene file: %s", e)
        return []


class PPI_Network:
    def __init__(self, result_dir):
        self.result_dir = result_dir
        os.makedirs(result_dir, exist_ok=True)

    def get_string_data(self, gene_list, taxon_id="9606", min_confidence=700):
        """Fetch STRING interactions with adjustable confidence cutoff."""
        logger = logging.getLogger()
        base_url = "https://string-db.org/api"
        format_type = "tsv"
        api_method = "network"
        api_params = {
            "identifiers": "%0d".join(gene_list),
            "species": taxon_id,
            "required_score": min_confidence,
            "caller_identity": "PPINetworkAnalysis"
        }
        full_url = f"{base_url}/{format_type}/{api_method}"
        try:
            api_response = requests.post(full_url, data=api_params, verify=False)
            api_response.raise_for_status()
            return api_response.text
        except requests.RequestException as e:
            logger.error("[ERROR]  Error fetching STRING interactions: %s", e)
            return

    def construct_network(self, gene_list, min_confidence=700):
        
        logger = logging.getLogger()
        logger.info("[INFO]  Building PPI network...")
        string_response = self.get_string_data(gene_list, min_confidence=min_confidence)
        Graph = nx.Graph()
        if string_response:
            try:
                interaction_df = pd.read_csv(StringIO(string_response), sep="\t")
                for _, edge_row in interaction_df.iterrows():
                    Graph.add_edge(edge_row["preferredName_A"], edge_row["preferredName_B"], weight=edge_row["score"])
                logger.info("[INFO]  PPI: %d nodes, %d edges", Graph.number_of_nodes(), Graph.number_of_edges())
            except pd.errors.EmptyDataError:
                logger.warning("[WARN]  No interactions from STRING API")

        graphml_file = os.path.join(self.result_dir, "ppi_network.graphml")
        try:
            nx.write_graphml(Graph, graphml_file)
            logger.info("[INFO]  GraphML: %s", graphml_file)
        except OSError as e:
            logger.error("[ERROR]  Error saving PPI network: %s", e)
        return Graph

    def get_hub_gene(self, ppi_graph):
        
        logger = logging.getLogger()
        if not ppi_graph.nodes():
            logger.warning("[WARN]  Network is empty, no hub genes")
            return pd.DataFrame(columns=["Gene", "Degree"])
        node_degrees = dict(ppi_graph.degree())
        mean_degree = sum(node_degrees.values()) / len(node_degrees) if node_degrees else 0
        high_degree_genes = [gene for gene, degree in node_degrees.items() if degree > mean_degree]
        logger.info("[INFO] ✔ %d hub genes (degree > mean)", len(high_degree_genes))
        hub_data = pd.DataFrame(
            [(gene, degree) for gene, degree in node_degrees.items() if gene in high_degree_genes],
            columns=["Gene", "Degree"]
        )
        hub_csv_file = os.path.join(self.result_dir, "hub_genes.csv")
        try:
            hub_data.to_csv(hub_csv_file, index=False)
            logger.info("[INFO]  Hub genes: %s", hub_csv_file)
        except OSError as e:
            logger.error("[ERROR]  Error saving hub genes: %s", e)
            return hub_data
        return hub_data

    def extract_hub_subgraph(self, ppi_graph, hub_genes):
        
        nodes_to_keep = set(hub_genes)
        for hub in hub_genes:
            nodes_to_keep.update(ppi_graph.neighbors(hub))
        return ppi_graph.subgraph(nodes_to_keep).copy()

    def render_network(self, ppi_graph, imp_genes, html_file, title="Protein–Protein Interaction (PPI) Network"):
        
        logger = logging.getLogger()
        if not ppi_graph.nodes():
            logger.warning("[WARN]  Empty network, skipping visualization")
            return

        logger.info("[INFO]  PPI network visualization...")

        
        vis_network = Network(
            height="700px", width="100%", bgcolor="#FFFFFF", font_color="black",
            notebook=True, cdn_resources='in_line'
        )

        
        vis_network.force_atlas_2based(
            gravity=-30,
            central_gravity=0.01,
            spring_length=100,
            spring_strength=0.05,
            damping=0.9,
            overlap=0
        )

        
        for node_id in ppi_graph.nodes():
            if node_id in imp_genes:  
                node_color = "#FF0000"  
            else:
                node_color = "#4682B4"  

            vis_network.add_node(
                node_id,
                shape="dot",
                label=node_id if node_id in imp_genes else "",
                color=node_color,
                size=20,
                title=f"{node_id}\nDegree: {ppi_graph.degree(node_id)}"
            )

        
        for source_node, target_node in ppi_graph.edges():
            vis_network.add_edge(source_node, target_node, color="#C0C0C0", width=1)

        if html_file is None:
            html_file = os.path.join(self.result_dir, "ppi_network_interactive.html")

        html_path = str(Path(html_file))
        try:
            vis_html = vis_network.generate_html(notebook=True)

            # Inject Title + Legend
            legend_html = f"""
            <div style="text-align:center; font-family:Arial; margin-bottom:10px;">
                <h2 style="color:#333;">{title}</h2>
                <div style="display:flex; justify-content:center; gap:30px; margin-top:5px;">
                    <div><span style="display:inline-block; width:15px; height:15px; background:#FF0000; border-radius:50%; margin-right:5px;"></span> Hub Gene</div>
                    <div><span style="display:inline-block; width:15px; height:15px; background:#4682B4; border-radius:50%; margin-right:5px;"></span> Other Gene</div>
                </div>
            </div>
            """
            vis_html = vis_html.replace("<body>", f"<body>{legend_html}", 1)

            with open(html_file, "w", encoding="utf-8") as file_handle:
                file_handle.write(vis_html)

            logger.info("[INFO]  HTML with title and legend: %s", html_file)
        except (OSError, UnicodeEncodeError) as vis_error:
            logger.error("[ERROR]  Error saving PPI network: %s", vis_error)
            return

        if boolean_track:
            try:
                iframe_html = f'<iframe src="{html_path}" width="100%" height="750px" frameborder="0"></iframe>'
                display(HTML(iframe_html))
                logger.info("[INFO] Presenting PPI network with title & legend")
            except Exception as e:
                logger.error("[ERROR]  Failed to display PPI network: %s", str(e))


def PPI_Analysis(gene_csv_file, min_conf=700, hub_only=True):
    logger = logging.getLogger()
    gene_name = gene_path(gene_csv_file)
    if not gene_name:
        logger.error("[ERROR]  No genes loaded")
        return
    res_dir = "output"
    ppi_builder = PPI_Network(res_dir)
    network_graph = ppi_builder.construct_network(gene_name, min_confidence=min_conf)
    if not network_graph.nodes():
        logger.error("[ERROR]  No network built")
        return
    hub_gene_df = ppi_builder.get_hub_gene(network_graph)
    hub_gene_name = hub_gene_df["Gene"].tolist()

    
    graph_to_render = network_graph
    if hub_only and hub_gene_name:
        graph_to_render = ppi_builder.extract_hub_subgraph(network_graph, hub_gene_name)

    ppi_builder.render_network(graph_to_render, hub_gene_name, None)