## Retrieve uniprot protein information ##
import requests

## Using gene names ##
def fetch_gene_info(gene_list):
    gene_annotations = {}
    base_url = "http://mygene.info/v3/query"
    for gene in gene_list:
        params = {
            "q": gene,
            "fields": "symbol, name, summary",
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            gene_data = response.json()
            if gene_data["hits"]:  # If gene found
                first_hit = gene_data["hits"][0]
                gene_symbol = first_hit.get("symbol", "N/A")
                gene_name = first_hit.get("name", "N/A")
                gene_summary = first_hit.get("summary", "N/A")

                gene_annotations[gene] = {"Gene": gene_symbol,
                                          "Gene_name": gene_name,
                                          "Gene_summary": gene_summary}
            else:
                print(f"No information found for gene: {gene}\n")
                gene_annotations[gene] = {"Gene": "Not found",
                                          "Gene_name": "Not found",
                                          "Gene_summary": "Not found"}
        else:
            print(f"Failed to fetch information for gene: {gene}\n")
            gene_annotations[gene] = "Search not made."
            
    return gene_annotations

## Using uniprot id (needs work) ##
def fetch_protein_info(protein_list):
    base_url = "https://www.uniprot.org/uniprotkb"
    for protein in protein_list:
        response = requests.get(f"{base_url}/{protein}")
        if response.status_code == 200:
            protein_data = response.json()
            if 'entry' in protein_data:
                entry = protein_data['entry']

                protein_id = entry.get('accession', 'N/A')
                protein_name = entry.get('protein', {}).get('recommendedName', {}).get('fullName', 'N/A')
                protein_function = entry.get('comments', [{'text': 'No function information available'}])[0]['text']

                print(f"Protein ID: {protein_id}")
                print(f"Protein Name: {protein_name}")
                print(f"Function: {protein_function}\n")
            else:
                print(f"No information found for protein: {protein}\n")
        else:
            print(f"Failed to fetch information for protein: {protein}\n")