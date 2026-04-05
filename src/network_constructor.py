import os
import networkx as nx
import pandas as pd
from pyvis.network import Network
from IPython.display import IFrame, display
import logging

def construct_circrna_mirna_mrna_network(results, strong_hits, all_mirnas, temp_dir, output_dir):
    logger = logging.getLogger()
    G = nx.DiGraph()
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    
    for circ in results.keys():
        match_file = os.path.join(temp_dir, f"{circ}_strong_medium_matches.csv")
        if os.path.exists(match_file):
            try:
                df = pd.read_csv(match_file)
                for _, row in df.iterrows():
                    c, m = row['circ_id'], row['mirna_id']
                    G.add_node(c, type='circRNA')
                    G.add_node(m, type='miRNA')
                    G.add_edge(c, m, interaction='circRNA→miRNA')
            except Exception as e:
                logger.error(f"[ERROR] Failed: {match_file} | {e}")

    
    overlap_file = os.path.join(temp_dir, "overlapping_mrnas.csv")
    if not os.path.exists(overlap_file):
        logger.warning("[WARN]  No miRNA–mRNA file found")
        return None

    try:
        overlap_df = pd.read_csv(overlap_file)
        for _, row in overlap_df.iterrows():
            m, g = row['mirna'], row['gene']
            G.add_node(m, type='miRNA')
            G.add_node(g, type='mRNA')
            G.add_edge(m, g, interaction='miRNA→mRNA')
    except Exception as e:
        logger.error("[ERROR]  Overlap load failed: %s", e)
        return None

    logger.info("[INFO]  Network: %d nodes, %d edges (%d circRNAs, %d miRNAs, %d mRNAs)",
                G.number_of_nodes(), G.number_of_edges(),
                len([n for n, d in G.nodes(data=True) if d.get('type') == 'circRNA']),
                len([n for n, d in G.nodes(data=True) if d.get('type') == 'miRNA']),
                len([n for n, d in G.nodes(data=True) if d.get('type') == 'mRNA']))


    graphml_path = os.path.join(output_dir, "circrna_mirna_mrna_network.graphml")
    try:
        nx.write_graphml(G, graphml_path)
        logger.info("[INFO]  GraphML: %s", graphml_path)
    except Exception as e:
        logger.error("[ERROR]  GraphML save failed: %s", e)

    
    
    
    net = Network(height="800px", width="100%", directed=True, notebook=True, cdn_resources="in_line")

    
    for node, attr in G.nodes(data=True):
        ntype = attr.get('type', 'unknown')

        if ntype == 'circRNA':
            color, shape, size, font_size, physics = 'lightpink', 'diamond', 45, 20, True
        elif ntype == 'miRNA':
            
            color, shape, size, font_size, physics = 'darkorange', 'triangle', 40, 18, False
        elif ntype == 'mRNA':
            color, shape, size, font_size, physics = 'lightgreen', 'ellipse', 25, 14, True
        else:
            color, shape, size, font_size, physics = 'lightgray', 'dot', 20, 12, True

        net.add_node(
            node,
            label=node,
            title=f"{ntype}: {node}",
            color=color,
            shape=shape,
            size=size,
            font={"size": font_size, "vadjust": 0},
            physics=physics
        )

    
    for src, tgt, attr in G.edges(data=True):
        inter = attr.get('interaction', '')

        if inter == 'circRNA→miRNA':
            color, width, dashed = 'darkblue', 4, True
        elif inter == 'miRNA→mRNA':
            color, width, dashed = 'gray', 1, False
        else:
            color, width, dashed = 'lightgray', 1, False

        net.add_edge(
            src, tgt,
            title=f"{src} → {tgt}",
            dashes=dashed,
            color=color,
            width=width
        )

    
    
    
    net.set_options("""{
      "physics": {
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 150,
          "springConstant": 0.05
        },
        "minVelocity": 0.75,
        "solver": "forceAtlas2Based"
      },
      "interaction": {
        "hover": true,
        "tooltipDelay": 0,
        "multiselect": true,
        "navigationButtons": false,
        "keyboard": true,
        "dragNodes": true
      },
      "nodes": {
        "font": {
          "size": 16,
          "strokeWidth": 0,
          "bold": "16px Arial black",
          "vadjust": 0
        }
      }
    }""")

    
    
    
    html_path = os.path.join(output_dir, "circrna_mirna_mrna_network.html")
    try:
        html_str = net.generate_html()

        
        header_html = """
        <div style='text-align:center; font-size:22px; font-weight:bold; padding:10px;'>
            circRNA–miRNA–mRNA Regulatory Network
        </div>

        <div style='text-align:center; font-size:16px; padding:5px;'>
            <span style="color:lightpink; font-weight:bold;">◆ circRNA</span> &nbsp;&nbsp;
            <span style="color:darkorange; font-weight:bold;">▲ miRNA</span> &nbsp;&nbsp;
            <span style="color:lightgreen; font-weight:bold;">● mRNA</span>
        </div>
        """
        html_str = header_html + html_str

        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_str)

        logger.info("[INFO]  HTML saved: %s", html_path)
    except Exception as e:
        logger.error("[ERROR]  HTML save failed: %s", e)
        return G

    
    try:
        display(IFrame(html_path, width=950, height=750))
    except Exception:
        pass

    return G