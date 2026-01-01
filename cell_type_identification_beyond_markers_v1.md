# Cell Type Identification Beyond Predefined Marker Lists

## The Challenge

**Problem**: Gene lists from studies may contain genes not explicitly listed in skill marker databases, yet these genes can still indicate specific cell types through:
1. Pathway/functional associations
2. Gene Ontology annotations
3. Literature-based evidence
4. Expression database patterns
5. Co-expression networks

**Solution**: Multi-method approach combining computational tools, databases, and reasoning strategies.

---

## CRITICAL: Special Considerations for DEG Analysis and TME Characterization

### The Top-N Gene Fallacy

**PROBLEM:** When analyzing large differential expression gene (DEG) datasets (e.g., 3,000+ genes), selecting only the top 20-50 genes by fold-change is INSUFFICIENT for comprehensive analysis, especially for tumor microenvironment (TME) characterization.

**Why Top-N Selection Fails:**

1. **Misses rare but critical cell types**
   - Example: Plasma cells may have moderate FC but are clinically significant
   - TLS (tertiary lymphoid structures) markers may not be in top 30

2. **Biases toward abundant cell populations**
   - Cancer cells (high abundance) dominate top genes
   - Rare immune populations (low abundance) are missed

3. **Loses pathway complexity**
   - May capture cancer pathways but miss immune checkpoints
   - May see T cell recruitment signals but miss T cell exclusion

4. **Misses bidirectional changes**
   - Top upregulated ≠ complete TME gained populations
   - Top downregulated ≠ complete TME lost populations

5. **Ignores functional modules**
   - May get 1 neutrophil marker but miss complete neutrophil exclusion
   - May see CXCL13 but miss entire TLS program

### Stratified Gene Selection Strategy for TME Analysis

**MANDATORY for comprehensive TME characterization:**

#### Upregulated Gene Selection (Gained Populations):

```python
def select_upregulated_tme_genes(deg_df, min_genes=100):
    """
    Stratified selection to capture ALL TME components in upregulated genes.
    """
    selected = []
    
    # 1. Top by fold-change (cancer cells, major changes)
    top_fc = deg_df[deg_df['log2FC'] > 1].nlargest(40, 'log2FC')
    selected.append(top_fc)
    
    # 2. Highly expressed (abundant populations regardless of FC)
    high_expr = deg_df[(deg_df['log2FC'] > 0.5) & 
                       (deg_df['aveEXP'] > 7)].nlargest(20, 'aveEXP')
    selected.append(high_expr)
    
    # 3. Immune markers (adaptive and innate)
    immune_patterns = ['CD', 'HLA', 'IL', 'CXCL', 'CCL', 'GZMB', 'PRF1', 
                      'IFNG', 'TNF', 'CTLA', 'PDCD', 'LAG3', 'TIGIT', 
                      'ICOS', 'CD28', 'FOXP3', 'TBET', 'GATA3', 'RORC']
    immune = deg_df[deg_df['name'].str.contains('|'.join(immune_patterns), 
                    case=False, na=False) & (deg_df['log2FC'] > 0)].head(30)
    selected.append(immune)
    
    # 4. Stromal/CAF markers
    stromal_patterns = ['COL', 'FN1', 'ACTA2', 'FAP', 'PDGFR', 'VIM', 
                       'POSTN', 'THY1', 'S100A4', 'TWIST', 'SNAI']
    stromal = deg_df[deg_df['name'].str.contains('|'.join(stromal_patterns), 
                     case=False, na=False) & (deg_df['log2FC'] > 0)].head(15)
    selected.append(stromal)
    
    # 5. ECM remodeling
    ecm_patterns = ['MMP', 'TIMP', 'LOXL', 'PLOD', 'ADAM', 'CTHRC', 'SPP1']
    ecm = deg_df[deg_df['name'].str.contains('|'.join(ecm_patterns), 
                 case=False, na=False) & (deg_df['log2FC'] > 0)].head(15)
    selected.append(ecm)
    
    # 6. Angiogenesis
    angio_patterns = ['VEGF', 'FLT', 'KDR', 'ANGPT', 'TEK', 'NRP', 'PLGF']
    angio = deg_df[deg_df['name'].str.contains('|'.join(angio_patterns), 
                   case=False, na=False) & (deg_df['log2FC'] > 0)].head(10)
    selected.append(angio)
    
    # 7. Hypoxia/metabolism
    hypoxia_patterns = ['CA9', 'BNIP3', 'LDHA', 'SLC2A1', 'ENO', 'PDK', 
                       'HIF', 'EGLN', 'PFKFB']
    hypoxia = deg_df[deg_df['name'].str.contains('|'.join(hypoxia_patterns), 
                     case=False, na=False) & (deg_df['log2FC'] > 0)].head(10)
    selected.append(hypoxia)
    
    # Combine and remove duplicates
    combined = pd.concat(selected).drop_duplicates(subset='name')
    return combined
```

#### Downregulated Gene Selection (Lost Populations):

```python
def select_downregulated_tme_genes(deg_df, min_genes=100):
    """
    Stratified selection to capture ALL lost cell types.
    """
    selected = []
    
    # 1. Top by fold-change (major losses)
    top_fc = deg_df[deg_df['log2FC'] < -1].nsmallest(40, 'log2FC')
    selected.append(top_fc)
    
    # 2. Highly expressed lost populations
    high_expr = deg_df[(deg_df['log2FC'] < -0.5) & 
                       (deg_df['aveEXP'] > 7)].nsmallest(20, 'log2FC')
    selected.append(high_expr)
    
    # 3. Immune markers lost (critical for TME!)
    immune_patterns = ['CD', 'HLA', 'IL', 'CXCL', 'CCL', 'FCGR', 'CSF', 
                      'ITGA', 'DEFA', 'S100A', 'MPO', 'ELANE', 'LYZ',
                      'GZMB', 'PRF1', 'NKG7', 'GNLY']
    immune = deg_df[deg_df['name'].str.contains('|'.join(immune_patterns), 
                    case=False, na=False) & (deg_df['log2FC'] < 0)].head(40)
    selected.append(immune)
    
    # 4. Endothelial markers
    endo_patterns = ['PECAM', 'VWF', 'CDH5', 'CLDN5', 'PLVAP', 'ESM1', 
                    'FLT1', 'TEK', 'ESAM', 'EMCN']
    endo = deg_df[deg_df['name'].str.contains('|'.join(endo_patterns), 
                  case=False, na=False) & (deg_df['log2FC'] < 0)].head(15)
    selected.append(endo)
    
    # 5. Muscle/pericytes
    muscle_patterns = ['ACTA2', 'MYH', 'TAGLN', 'DES', 'CNN', 'MYLK', 
                      'RGS5', 'PDGFR', 'MCAM', 'NOTCH3']
    muscle = deg_df[deg_df['name'].str.contains('|'.join(muscle_patterns), 
                    case=False, na=False) & (deg_df['log2FC'] < 0)].head(15)
    selected.append(muscle)
    
    # 6. Normal epithelial (if cancer sample)
    epithelial_patterns = ['KRT', 'CDH1', 'EPCAM', 'CLDN', 'TJP', 'OCLN']
    epithelial = deg_df[deg_df['name'].str.contains('|'.join(epithelial_patterns), 
                        case=False, na=False) & (deg_df['log2FC'] < 0)].head(10)
    selected.append(epithelial)
    
    # Combine and remove duplicates
    combined = pd.concat(selected).drop_duplicates(subset='name')
    return combined
```

### Recommended Gene Numbers for Different Analyses

| Analysis Type | Minimum Genes | Recommended | Rationale |
|--------------|---------------|-------------|-----------|
| **Simple cell type ID** | 20-30 | 50 | Basic identification |
| **TME characterization** | 100 | 200-300 | Capture all compartments |
| **Immune landscape** | 50 | 100-150 | All immune populations |
| **Cancer + TME** | 150 | 250-400 | Cancer + all TME components |
| **Pathway analysis only** | 30 | 50-100 | Functional modules |

### Critical TME Components to Capture

**ALWAYS check for these in TME analysis:**

#### 1. Adaptive Immune Cells:
- **T cells:** CD3D/E/G, CD4, CD8A/B, FOXP3, GATA3, TBX21, RORC
- **B cells:** CD19, CD79A/B, MS4A1, PAX5
- **Plasma cells:** SDC1, MZB1, JCHAIN, IGHG1/2/3/4
- **Tertiary lymphoid structures (TLS):** CXCL13, CCL19, CCL21, LAMP3

#### 2. Innate Immune Cells:
- **Neutrophils:** MPO, ELANE, DEFA1/3, S100A8/9/12, FCGR3B
- **Monocytes:** CD14, FCGR1A, S100A8/9, LYZ, VCAN
- **Macrophages:** CD68, CD163, MSR1, MRC1, MARCO, C1QA/B/C
- **Dendritic cells:** CD1C, CLEC9A, CLEC10A, FLT3, ITGAX, LAMP3
- **NK cells:** NCAM1, FCGR3A, KLRD1, KLRB1, GNLY, NKG7
- **Mast cells:** TPSAB1, TPSB2, CPA3, KIT, HDC
- **Eosinophils:** RNASE2/3, EPX, PRG2

#### 3. Immune Checkpoints:
- **Inhibitory:** PDCD1, CTLA4, LAG3, TIGIT, HAVCR2, BTLA, ADORA2A
- **Stimulatory:** CD28, ICOS, CD27, CD40LG, TNFRSF4, TNFRSF9

#### 4. Cytokines/Chemokines:
- **Pro-inflammatory:** IL1B, IL6, TNF, IFNG, IL12A/B, IL18
- **Anti-inflammatory:** IL10, TGFB1, IL4, IL13
- **Chemokines:** CXCL8/9/10/11/13, CCL2/3/4/5/18/19/20/22

#### 5. Antigen Presentation:
- **MHC-I:** HLA-A/B/C, B2M
- **MHC-II:** HLA-DRA/DRB1, HLA-DQA1/DQB1, HLA-DPA1/DPB1, CD74

#### 6. Stromal Components:
- **CAFs:** FAP, ACTA2, COL1A1/A2, FN1, PDGFRA/B, POSTN, THY1
- **ECM:** Collagens (COL1-28), MMPs, TIMPs, LOXL1/2
- **Endothelial:** PECAM1, VWF, CDH5, CLDN5, FLT1, KDR
- **Pericytes:** RGS5, PDGFRB, MCAM, NOTCH3, DES

#### 7. Functional Modules:
- **Angiogenesis:** VEGFA/B/C, FLT1, KDR, ANGPT1/2, TEK
- **Hypoxia:** CA9, BNIP3, HIF1A, EGLN1/3, LDHA, SLC2A1
- **EMT:** TWIST1/2, SNAI1/2, ZEB1/2, VIM, CDH1/2
- **Proliferation:** MKI67, TOP2A, PCNA, CDK1, CCNB1

### Real-World Example: Lung Adenocarcinoma TME

**Scenario:** 3,167 DEGs from tumor vs normal lung tissue

**Top 30 genes approach (WRONG):**
- Captured: Cancer cells (CEACAM5, SPP1, MUC5B)
- **MISSED:** Plasma cells, TLS, Tregs, T cell exclusion, neutrophil exclusion, macrophage loss

**Stratified 224 genes approach (CORRECT):**
- Captured: Cancer cells + Plasma cells + TLS + Tregs + T cell exclusion + Complete innate immune loss
- **Result:** Identified "Immune Paradox" phenotype (TLS+ but T-cell-excluded)
- **Clinical impact:** Changed recommendation from PD-1 monotherapy (won't work) to combination vascular normalization + immunotherapy

### TME Phenotype Classification Framework

Based on comprehensive gene selection, classify TME into:

#### 1. "Hot" / Inflamed TME
**Markers present:**
- ✅ CD8+ T cells (CD8A, GZMB, PRF1)
- ✅ IFN-γ signature (CXCL9, CXCL10, IFNG)
- ✅ M1 macrophages (NOS2, IL1B, TNF)
- ✅ Activated DCs (LAMP3, CCR7)
- ✅ High MHC-I/II

**Treatment:** PD-1/PD-L1 inhibitors likely effective

#### 2. "Cold" / Excluded TME
**Markers present:**
- ❌ Few/no CD8+ T cells
- ❌ Low chemokines
- ❌ Few immune cells
- ❌ CAFs/stroma high (physical barrier)

**Treatment:** Needs combination therapy (vascular normalization, radiation)

#### 3. "Immune Paradox" / Marginalized TME (LITERATURE-VALIDATED!)
**Markers present:**
- ✅ TLS present (CXCL13, plasma cells: SDC1, MZB1, JCHAIN)
- ✅ T cell recruitment signals (CXCL9/10)
- ❌ But CD8+ T cells EXCLUDED (PRF1-, NKG7-, GNLY-)
- ⚠️ Tregs present (FOXP3+)
- ⚠️ Treg recruitment (CCL22+)
- ❌ Innate immunity absent (neutrophils, macrophages, DCs all lost)
- ❌ Vascular dysfunction (CLDN5 loss, endothelial markers lost)

**Literature Support:**
- **Joshi et al., Immunity 2015:** Tregs in TLS suppress anti-tumor T cells in lung adenocarcinoma mouse model; Treg depletion rescues immunity
- **MedRxiv 2020:** TLS+ with plasma cells but CD8+ excluded in early-stage lung adenocarcinoma; follicular Tregs in TLS explain exclusion
- **J Hematol Oncol 2025:** Immature TLS (PFL-TLS) with Tregs are immunosuppressive and predict poor ICI response
- **Mechanism:** Tregs suppress T cell priming in TLS + vascular abnormalities block T cell entry + immature TLS lack functional germinal centers

**Treatment (Evidence-Based):** 
1. **Vascular normalization:** Bevacizumab (anti-VEGF) - enable T cell entry, restore CLDN5
2. **Treg depletion:** Anti-CCR4 (mogamulizumab) or metronomic cyclophosphamide
3. **Immunotherapy + chemotherapy:** Atezolizumab + carboplatin/pemetrexed + bevacizumab (FDA-approved: IMpower150 regimen)
4. **TLS maturation (experimental):** STING agonist (intratumoral) - activate DCs, boost TLS maturation, recruit innate immunity

**Key Insight:** TLS presence ≠ functional anti-tumor immunity. Must assess TLS maturity, composition (Treg:CD8 ratio), and vascular competence. PD-1 monotherapy will fail because T cells are EXCLUDED, not exhausted.

#### 4. "Suppressed" / Exhausted TME
**Markers present:**
- ✅ CD8+ T cells PRESENT
- ✅ But exhaustion markers HIGH (PDCD1, CTLA4, LAG3, TIGIT)
- ⚠️ M2 macrophages (CD163, MRC1)
- ⚠️ MDSCs
- ⚠️ Tregs

**Treatment:** PD-1/PD-L1 + CTLA-4 combination

### Validation Checklist for TME Analysis

Before finalizing TME characterization, verify:

- [ ] **Coverage:** All major immune populations checked (T, B, NK, neutrophils, macrophages, DCs)
- [ ] **Bidirectional:** Both upregulated AND downregulated genes analyzed
- [ ] **Rare cells:** Plasma cells, Tregs, DCs explicitly checked (may have moderate FC)
- [ ] **Functional modules:** Checkpoints, cytokines, chemokines, MHC assessed
- [ ] **Phenotype coherence:** Do findings fit a known TME phenotype?
- [ ] **Clinical actionability:** Does phenotype guide treatment selection?

### Common Pitfalls to Avoid

**❌ DON'T:**
1. Use only top 20-50 genes for TME analysis
2. Ignore downregulated genes (lost populations are critical!)
3. Assume high FC = high importance (rare cells may have moderate FC)
4. Miss rare but critical cell types (plasma cells, Tregs, DCs)
5. Classify as simply "hot" or "cold" (TME is often nuanced)
6. Ignore spatial patterns (center vs edge, stroma vs tumor)

**✅ DO:**
1. Use stratified sampling (100-300 genes minimum)
2. Analyze BOTH upregulated AND downregulated
3. Explicitly check for rare cell types
4. Check for paradoxes (recruitment signals but cells absent)
5. Use nuanced phenotype classification
6. Consider spatial heterogeneity when interpreting

---

## Method 1: Enrichment-Based Cell Type Inference

### Concept
Instead of relying on individual marker genes, use **gene set enrichment** against cell-type-specific gene expression databases.

### Databases for Cell Type Enrichment

**1. Human Gene Atlas / Mouse Gene Atlas**
- Tissue/cell line expression profiles
- Available via Enrichr, GSEA
- Example: Gene list enriched in "CD8+ T cells" or "Fibroblasts"

**2. ARCHS4 Tissues / ARCHS4 Cell Lines**
- Massive expression compendium
- Cell-type and tissue-specific signatures
- Available via Enrichr

**3. Tabula Sapiens / Tabula Muris**
- Single-cell RNA-seq atlases
- Cell-type-specific gene sets from healthy tissues
- Comprehensive cell type coverage

**4. Human Protein Atlas**
- Cell-type-specific protein expression
- RNA expression across cell types
- Can query genes individually

**5. PanglaoDB**
- Single-cell marker database
- Cell type marker genes across species
- Web interface and downloadable

**6. CellMarker 2.0**
- Curated cell markers from literature
- Tissue and cell type specific
- Covers normal and cancer cells

**7. Azimuth Reference Atlases**
- PBMC reference (immune cells)
- Lung, kidney, motor cortex, etc.
- Can map gene signatures to reference

### Implementation Strategy

```python
def infer_cell_type_by_enrichment(gene_list):
    """
    Use enrichment against cell type databases to infer cell identity.
    """
    import gseapy as gp
    
    # Multiple cell type databases
    cell_type_libraries = [
        'Human_Gene_Atlas',           # Broad tissue/cell types
        'ARCHS4_Tissues',             # Expression-based
        'Tabula_Sapiens',             # Single-cell based
        'ARCHS4_Cell-lines',          # Cell line expression
        'Jensen_TISSUES'              # Text-mining based
    ]
    
    results = {}
    for library in cell_type_libraries:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=library,
                organism='human',
                outdir=None
            )
            results[library] = enr.results.head(5)
        except:
            continue
    
    # Consensus across databases
    return identify_consensus_cell_type(results)
```

**Advantages:**
- Works with ANY genes, not just known markers
- Statistical significance (p-values)
- Considers entire gene expression pattern
- Multiple databases provide validation

**Example:**
```
Gene list: [CTHRC1, WNT5A, IL13]
↓
Enrichr → Human_Gene_Atlas
↓
Results:
1. Fibroblasts (p = 1e-3) - CTHRC1 enriched
2. Macrophages (p = 5e-3) - WNT5A enriched  
3. Th2 cells (p = 1e-2) - IL13 enriched
↓
Interpretation: Mixed TME signature
```

---

## Method 2: Gene Ontology (GO) and Functional Annotation

### Concept
Use GO terms, KEGG pathways, and functional annotations to infer biological processes, then map processes to cell types.

### Strategy

**Step 1: Functional Enrichment**
```python
def functional_to_cell_type(gene_list):
    """
    Map GO terms and pathways to cell types.
    """
    # Get GO enrichment
    go_results = run_go_enrichment(gene_list, ['GO_Biological_Process_2023'])
    
    # Map GO terms to cell types
    go_to_cell_type = {
        'T cell activation': 'T cells',
        'B cell proliferation': 'B cells',
        'natural killer cell mediated cytotoxicity': 'NK cells',
        'macrophage activation': 'Macrophages',
        'collagen fibril organization': 'Fibroblasts',
        'angiogenesis': 'Endothelial cells',
        'antigen processing and presentation': 'APCs (DCs, macrophages)',
        'mast cell activation': 'Mast cells',
        'neutrophil chemotaxis': 'Neutrophils',
        'type 2 immune response': 'Th2 cells, eosinophils',
        # ... many more mappings
    }
    
    # Infer cell types from enriched processes
    cell_types = []
    for go_term in go_results['Term']:
        for pattern, cell_type in go_to_cell_type.items():
            if pattern.lower() in go_term.lower():
                cell_types.append(cell_type)
    
    return cell_types
```

**Example Mappings:**

| GO Term / Process | Likely Cell Type |
|-------------------|------------------|
| T cell activation | T cells |
| Regulation of cytokine production | Immune cells (broad) |
| Extracellular matrix organization | Fibroblasts, CAFs |
| Angiogenesis | Endothelial cells |
| Phagocytosis | Macrophages, neutrophils |
| Antigen presentation | DCs, macrophages, B cells |
| NK cell cytotoxicity | NK cells |
| Immunoglobulin production | Plasma cells |
| Mast cell degranulation | Mast cells |
| Epithelial cell proliferation | Epithelial cells |
| Wound healing | Fibroblasts, macrophages |

**Step 2: Pathway-to-Cell-Type Mapping**

| KEGG Pathway | Likely Cell Type |
|--------------|------------------|
| T cell receptor signaling | T cells |
| B cell receptor signaling | B cells |
| NK cell mediated cytotoxicity | NK cells |
| Fc gamma R-mediated phagocytosis | Macrophages, neutrophils |
| Antigen processing and presentation | APCs |
| ECM-receptor interaction | Fibroblasts, CAFs |
| Focal adhesion | Fibroblasts, CAFs |
| VEGF signaling pathway | Endothelial cells |
| Tight junction | Epithelial cells |

---

## Method 3: Literature Mining and Database Queries

### Concept
For genes not in predefined lists, query literature and expression databases to find cell-type associations.

### Tools and Approaches

**1. PubMed / PubTator Central**
- Search: "[GENE] AND cell type"
- Search: "[GENE] AND expression"
- Look for phrases like "specifically expressed in" or "marker of"

**2. GeneCards**
- Comprehensive gene information
- Expression section shows tissue/cell type
- Links to Human Protein Atlas

**3. Human Protein Atlas**
- RNA expression by cell type (single-cell data)
- Protein expression by tissue
- Web interface: https://www.proteinatlas.org/
- API available for batch queries

**4. GTEx Portal**
- Tissue-specific expression
- Can infer cell type from tissue enrichment
- Web: https://gtexportal.org/

**5. Expression Atlas (EBI)**
- Baseline and differential expression
- Cell type and tissue specific
- https://www.ebi.ac.uk/gxa/

### Implementation

```python
def query_expression_databases(gene):
    """
    Query multiple databases for cell-type expression.
    """
    results = {}
    
    # Human Protein Atlas (if API available)
    hpa_data = query_hpa(gene)
    results['HPA'] = {
        'cell_types': extract_cell_types_from_hpa(hpa_data),
        'tissues': extract_tissues_from_hpa(hpa_data)
    }
    
    # GTEx
    gtex_data = query_gtex(gene)
    results['GTEx'] = {
        'top_tissues': get_top_expressed_tissues(gtex_data)
    }
    
    # Expression Atlas
    atlas_data = query_expression_atlas(gene)
    results['Expression_Atlas'] = {
        'cell_types': extract_cell_types(atlas_data)
    }
    
    return results

def infer_from_literature(gene):
    """
    Mine PubMed abstracts for cell type associations.
    """
    # Search PubMed
    query = f"{gene} AND (cell type OR expression OR marker)"
    abstracts = search_pubmed(query, max_results=20)
    
    # Extract cell type mentions
    cell_type_keywords = [
        'T cell', 'B cell', 'macrophage', 'fibroblast', 
        'NK cell', 'dendritic cell', 'neutrophil',
        'endothelial', 'epithelial', 'cancer-associated fibroblast'
    ]
    
    mentions = {}
    for abstract in abstracts:
        for cell_type in cell_type_keywords:
            if cell_type.lower() in abstract.lower():
                mentions[cell_type] = mentions.get(cell_type, 0) + 1
    
    # Rank by frequency
    return sorted(mentions.items(), key=lambda x: x[1], reverse=True)
```

---

## Method 4: Gene Co-Expression Networks

### Concept
Genes expressed together in the same cell type will co-express across samples. Use co-expression to infer cell type.

### Approach

**1. STRING Database**
- Protein-protein interaction and co-expression
- Query gene → get co-expressed genes
- Check if co-expressed genes are known markers

**2. COXPRESdb**
- Co-expression database
- Find genes co-expressed with your query gene
- Infer cell type from co-expressed markers

**3. ARCHS4 / SEEK**
- Co-expression from massive datasets
- Find similarly expressed genes
- Infer cell type from pattern

### Implementation

```python
def infer_from_coexpression(gene, known_markers):
    """
    Find cell type by identifying co-expressed known markers.
    """
    # Get top co-expressed genes (from STRING, COXPRESdb, etc.)
    coexpressed_genes = get_coexpressed_genes(gene, top_n=50)
    
    # Check overlap with known markers
    cell_type_scores = {}
    for cell_type, markers in known_markers.items():
        overlap = set(coexpressed_genes) & set(markers)
        if len(overlap) > 0:
            cell_type_scores[cell_type] = len(overlap)
    
    # Rank by overlap
    return sorted(cell_type_scores.items(), key=lambda x: x[1], reverse=True)
```

**Example:**
```
Query: CTHRC1
↓
STRING co-expressed: POSTN, FN1, COL1A1, FAP, ACTA2
↓
Check against known markers:
- CAF markers: ACTA2, FAP, COL1A1 ✓✓✓
- Endothelial: VWF, PECAM1 ✗
↓
Inference: CTHRC1 likely CAF marker
```

---

## Method 5: Single-Cell RNA-seq Reference Mapping

### Concept
Use single-cell reference atlases to map gene expression patterns to cell types.

### Tools

**1. Azimuth**
- Web-based: https://azimuth.hubmapconsortium.org/
- References: PBMC, lung, kidney, motor cortex, etc.
- Input: Gene expression → Output: Cell type

**2. CellTypist**
- Python package
- Pre-trained models for multiple tissues
- Fast automated annotation

```python
import celltypist

# Load model
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')

# Predict (requires expression matrix with genes as index)
predictions = celltypist.annotate(
    adata,  # AnnData object with your genes
    model=model,
    majority_voting=True
)
```

**3. SingleR (R package)**
- Reference-based annotation
- Uses correlation with reference profiles
- Multiple reference datasets available

```r
library(SingleR)

# Predict cell types
predictions <- SingleR(
    test = your_expression_data,
    ref = reference_data,
    labels = reference_data$cell_type
)
```

---

## Method 6: AI/LLM-Based Inference (When All Else Fails)

### Concept
Use Claude or other AI to reason about genes based on:
1. Gene name/symbol clues
2. Protein family membership
3. Known pathways
4. Biological knowledge

### Strategy

```
Prompt Template:

"Given the following genes: [GENE_LIST]

For each gene, determine:
1. What is the general function of this gene/protein?
2. Which cell types are likely to express this gene?
3. Is this gene a known marker for any cell type?
4. What biological processes involve this gene?

Use your knowledge of:
- Gene families (e.g., CD antigens, integrins, chemokines)
- Protein functions (e.g., transcription factors, surface receptors)
- Pathway participation
- Literature associations

Provide cell type inference with confidence level (high/medium/low)."
```

### Gene Name Clues

**Surface markers (CD antigens):**
- CD8A/CD8B → CD8+ T cells (high confidence)
- CD19, CD79A → B cells (high confidence)
- CD68, CD163 → Macrophages (high confidence)
- CD31 (PECAM1) → Endothelial (high confidence)

**Transcription factors:**
- FOXP3 → Regulatory T cells
- GATA3 → Th2 cells
- TBX21 (T-BET) → Th1 cells
- PAX5 → B cells
- PU.1 (SPI1) → Myeloid cells

**Chemokines/Cytokines:**
- CXCL9, CXCL10 → IFN-γ activated cells
- CCL18 → M2 macrophages
- IL13 → Th2 cells
- IFNG → Th1 cells, NK cells

**Structural proteins:**
- COL1A1, FN1 → Fibroblasts
- KRT (keratins) → Epithelial cells
- VIM (vimentin) → Mesenchymal cells

**Function-specific:**
- GZMB, PRF1 → Cytotoxic cells (CTL, NK)
- IGHG, JCHAIN → Plasma cells
- MKI67, TOP2A → Proliferating cells (any type)

---

## Integrated Multi-Method Workflow

### Step-by-Step Process

```python
def identify_cell_type_comprehensive(gene_list):
    """
    Comprehensive cell type identification using multiple methods.
    """
    results = {
        'methods': {},
        'consensus': None,
        'confidence': None
    }
    
    # Method 1: Check against known markers (fast)
    marker_result = check_known_markers(gene_list)
    results['methods']['markers'] = marker_result
    
    # Method 2: Enrichment analysis (high confidence)
    enrichment_result = enrichment_based_inference(gene_list)
    results['methods']['enrichment'] = enrichment_result
    
    # Method 3: GO/Pathway functional annotation
    functional_result = functional_annotation_inference(gene_list)
    results['methods']['functional'] = functional_result
    
    # Method 4: Expression database queries (for unknown genes)
    unknown_genes = identify_unknown_genes(gene_list, marker_result)
    if unknown_genes:
        expression_result = query_expression_databases_bulk(unknown_genes)
        results['methods']['expression_db'] = expression_result
    
    # Method 5: Co-expression analysis (if needed)
    if len(unknown_genes) > 0:
        coexpr_result = coexpression_inference(unknown_genes)
        results['methods']['coexpression'] = coexpr_result
    
    # Consensus building
    results['consensus'] = build_consensus(results['methods'])
    results['confidence'] = assess_confidence(results['methods'])
    
    return results


def build_consensus(method_results):
    """
    Build consensus cell type from multiple methods.
    """
    # Vote-based approach
    cell_type_votes = {}
    
    for method, result in method_results.items():
        if result and 'cell_type' in result:
            cell_types = result['cell_type']
            if isinstance(cell_types, str):
                cell_types = [cell_types]
            
            for ct in cell_types:
                cell_type_votes[ct] = cell_type_votes.get(ct, 0) + 1
    
    # Rank by votes
    ranked = sorted(cell_type_votes.items(), key=lambda x: x[1], reverse=True)
    
    # Return top consensus
    if ranked:
        return {
            'primary': ranked[0][0],
            'alternatives': [ct for ct, votes in ranked[1:3]],
            'votes': dict(ranked)
        }
    return None


def assess_confidence(method_results):
    """
    Assess confidence in cell type identification.
    """
    num_methods = len(method_results)
    num_concordant = count_concordant_methods(method_results)
    
    if num_concordant >= 3:
        return 'HIGH'
    elif num_concordant == 2:
        return 'MEDIUM'
    elif num_concordant == 1:
        return 'LOW'
    else:
        return 'UNCERTAIN'
```

---

## CRITICAL: Multi-Tissue Master Regulator Disambiguation

### The Problem

**Some master transcription factors create THE SAME transcriptional program in MULTIPLE tissues.**

This can cause framework to correctly identify the cellular program but incorrectly assign the tissue.

### Case Study: FOXI1 Paradox (Lesson Learned: Dec 31, 2025)

**FOXI1 drives identical programs in:**

| Tissue | Cell Type | Shared Markers | Unique Markers |
|--------|-----------|----------------|----------------|
| **Kidney** | Intercalated cells | FOXI1, ATP6V1B1, ATP6V1G3, BSND | SLC4A1/AE1, SLC26A4, CAII |
| **Lung** | Pulmonary ionocytes | FOXI1, ATP6V1B1, ATP6V1G3, BSND | **CFTR (extremely high!)** |
| **Inner Ear** | FORE cells | FOXI1, ATP6V1B1, ATP6V1G3, BSND | ATP6V0A4 |

**Real Example:**
- Gene list: `FOXI1, ATP6V1B1, ATP6V1G3, BSND, TFCP2L1, AQP6`
- Framework incorrectly identified: Kidney intercalated cells
- Correct answer: **Pulmonary ionocytes** (airway epithelium)
- Why wrong: Did not search multiple FOXI1 tissues, missed that CFTR discriminates

### Multi-Tissue Master Regulators - Watch List

**High-risk transcription factors that require tissue disambiguation:**

#### FOXI1
- **Tissues:** Kidney (ICs), Lung (ionocytes), Inner ear (endolymph), Epididymis
- **Discriminators:** 
  - CFTR (very high) → Lung ionocytes
  - SLC4A1/AE1 OR SLC26A4 → Kidney ICs
  - ATP6V0A4 → Inner ear

#### SOX2
- **Tissues:** Neural stem cells, Foregut endoderm, Tracheal epithelium, Retinal progenitors
- **Discriminators:** Context-specific co-factors (SOX1, SOX3 for neural; CDX2 for gut)

#### PAX6
- **Tissues:** Eye (lens, cornea, retina), Brain (neural progenitors), Pancreas (alpha cells)
- **Discriminators:**
  - SIX3, SIX6 → Eye development
  - PDX1, NKX6.1 → Pancreas
  - NEUROD1 → Retinal neurons

#### GATA3
- **Tissues:** T cells (Th2), Mammary gland, Kidney (nephric duct), Inner ear
- **Discriminators:**
  - CD3E, CD4, IL4 → T cells
  - KRT14, KRT18 → Mammary epithelium
  - LHX1, PAX2 → Kidney

#### FOXP3
- **Tissues:** Regulatory T cells, Some epithelial (controversial)
- **Discriminators:**
  - CD4, CD25, IL2RA, CTLA4 → Tregs
  - Without immune markers → Verify with literature

### Tissue Disambiguation Protocol

**MANDATORY when master TF with multi-tissue expression is identified:**

#### Step 1: Recognize Multi-Tissue TF

Check if identified transcription factor is known to operate in multiple tissues.

```python
MULTI_TISSUE_TFS = {
    'FOXI1': ['kidney', 'lung', 'inner_ear', 'epididymis'],
    'SOX2': ['neural', 'foregut', 'trachea', 'retina'],
    'PAX6': ['eye', 'brain', 'pancreas'],
    'GATA3': ['T_cells', 'mammary', 'kidney', 'inner_ear'],
    'FOXP3': ['T_cells', 'epithelial'],
    # Add more as discovered
}

def check_multitissue_tf(gene_list):
    """Check if gene list contains multi-tissue master TFs."""
    found_multitissue_tfs = []
    for gene in gene_list:
        if gene in MULTI_TISSUE_TFS:
            found_multitissue_tfs.append({
                'tf': gene,
                'tissues': MULTI_TISSUE_TFS[gene]
            })
    return found_multitissue_tfs
```

#### Step 2: Multi-Tissue Literature Search

**DO NOT only search one tissue!** Search ALL tissues where TF is expressed.

```python
def multitissue_literature_search(tf_name, core_genes):
    """
    Search literature for TF across all known tissues.
    """
    tissues = MULTI_TISSUE_TFS.get(tf_name, [])
    
    searches = []
    
    # Search each tissue specifically
    for tissue in tissues:
        query = f"{tf_name} {' '.join(core_genes[:3])} {tissue}"
        searches.append(query)
    
    # Broader searches
    searches.extend([
        f"{tf_name} all tissues expression",
        f"{tf_name} multiple cell types",
        f"{tf_name} {' '.join(core_genes[:5])} cell type"
    ])
    
    return searches
```

**Example for FOXI1 gene list:**
```
✅ Search 1: "FOXI1 ATP6V1B1 BSND kidney intercalated cells"
✅ Search 2: "FOXI1 ATP6V1B1 BSND lung airway ionocyte"
✅ Search 3: "FOXI1 ATP6V1B1 BSND inner ear"
✅ Search 4: "FOXI1 all tissues expression sites"
✅ Search 5: "FOXI1 cell types multiple"
```

#### Step 3: Build Tissue Discriminator Matrix

For each candidate tissue, identify which genes are **REQUIRED** vs **EXCLUDED**.

```python
TISSUE_DISCRIMINATORS = {
    'FOXI1': {
        'kidney_IC': {
            'required_any': ['SLC4A1', 'AE1', 'SLC26A4', 'CAII'],
            'excluded': ['CFTR_high'],  # CFTR rare/low in kidney
            'shared': ['ATP6V1B1', 'ATP6V1G3', 'BSND']
        },
        'lung_ionocyte': {
            'required_any': ['CFTR'],  # Extremely high!
            'excluded': ['SLC4A1', 'SLC26A4'],
            'shared': ['ATP6V1B1', 'ATP6V1G3', 'BSND']
        },
        'inner_ear': {
            'required_any': ['ATP6V0A4', 'SLC26A4'],
            'excluded': ['CFTR_high'],
            'shared': ['ATP6V1B1', 'ATP6V1G3', 'BSND']
        }
    }
}

def check_tissue_discriminators(gene_list, tf_name):
    """
    Check which tissue is most consistent with gene list.
    """
    if tf_name not in TISSUE_DISCRIMINATORS:
        return None
    
    tissue_scores = {}
    discriminators = TISSUE_DISCRIMINATORS[tf_name]
    
    for tissue, markers in discriminators.items():
        score = 0
        evidence = {'for': [], 'against': []}
        
        # Check required genes
        required = markers.get('required_any', [])
        has_required = any(gene in gene_list for gene in required)
        if has_required:
            score += 2
            evidence['for'].append(f"Has discriminator: {required}")
        
        # Check excluded genes (should NOT be present)
        excluded = markers.get('excluded', [])
        has_excluded = any(gene in gene_list for gene in excluded)
        if has_excluded:
            score -= 2
            evidence['against'].append(f"Has excluded marker: {excluded}")
        
        # Check shared genes
        shared = markers.get('shared', [])
        num_shared = sum(1 for gene in shared if gene in gene_list)
        score += num_shared * 0.5
        evidence['for'].append(f"Shared markers: {num_shared}/{len(shared)}")
        
        tissue_scores[tissue] = {
            'score': score,
            'evidence': evidence
        }
    
    return tissue_scores
```

#### Step 4: Missing Marker Analysis

**Check what's MISSING that should be present** for each candidate tissue.

```python
def analyze_missing_markers(gene_list, tf_name):
    """
    For each candidate tissue, check EXPECTED but MISSING markers.
    """
    if tf_name not in TISSUE_DISCRIMINATORS:
        return None
    
    missing_analysis = {}
    discriminators = TISSUE_DISCRIMINATORS[tf_name]
    
    for tissue, markers in discriminators.items():
        expected = markers.get('required_any', []) + markers.get('shared', [])
        missing = [gene for gene in expected if gene not in gene_list]
        
        missing_analysis[tissue] = {
            'expected_genes': expected,
            'missing_genes': missing,
            'percent_missing': len(missing) / len(expected) if expected else 0
        }
    
    return missing_analysis
```

**Example - FOXI1 gene list analysis:**

Gene list: `FOXI1, ATP6V1B1, ATP6V1G3, BSND, TFCP2L1`

```
Kidney IC:
  ✅ Has: ATP6V1B1, ATP6V1G3, BSND (shared)
  ❌ Missing: SLC4A1/AE1, SLC26A4, CAII (key discriminators!)
  📊 Missing: 3/7 expected genes (43%)

Lung Ionocyte:
  ✅ Has: ATP6V1B1, ATP6V1G3, BSND (shared)
  ❓ Missing: CFTR (but could be highest expresser, just not in list)
  📊 Missing: 1/4 expected genes (25%)

Conclusion: Lung ionocyte more likely (fewer missing expected markers)
```

#### Step 5: Build Tissue Confidence Matrix

Combine all evidence to assign confidence to each tissue.

```python
def build_tissue_confidence_matrix(
    gene_list, 
    tf_name, 
    literature_results,
    discriminator_scores,
    missing_analysis
):
    """
    Build comprehensive confidence assessment for each tissue.
    """
    tissues = MULTI_TISSUE_TFS.get(tf_name, [])
    confidence_matrix = {}
    
    for tissue in tissues:
        # Calculate composite score
        lit_score = literature_results.get(tissue, {}).get('relevance', 0)
        disc_score = discriminator_scores.get(tissue, {}).get('score', 0)
        missing_penalty = missing_analysis.get(tissue, {}).get('percent_missing', 0)
        
        composite = lit_score + disc_score - (missing_penalty * 2)
        
        # Assign confidence level
        if composite >= 3:
            confidence = 'HIGH'
        elif composite >= 1.5:
            confidence = 'MEDIUM-HIGH'
        elif composite >= 0:
            confidence = 'MEDIUM'
        else:
            confidence = 'LOW'
        
        confidence_matrix[tissue] = {
            'confidence': confidence,
            'composite_score': composite,
            'literature_score': lit_score,
            'discriminator_score': disc_score,
            'missing_penalty': missing_penalty
        }
    
    return confidence_matrix
```

#### Step 6: Report with Tissue Alternatives

**DO NOT report single tissue with VERY HIGH confidence if discriminators are missing!**

```python
def report_multitissue_results(
    tf_program,
    tissue_confidence_matrix,
    discriminators_present,
    discriminators_needed
):
    """
    Report cell type with appropriate tissue confidence.
    """
    # Sort tissues by confidence
    ranked_tissues = sorted(
        tissue_confidence_matrix.items(),
        key=lambda x: x[1]['composite_score'],
        reverse=True
    )
    
    top_tissue = ranked_tissues[0]
    alternatives = ranked_tissues[1:3]
    
    report = {
        'tf_program': tf_program,
        'program_confidence': 'VERY HIGH',  # We identified the TF program correctly
        'top_tissue': top_tissue[0],
        'tissue_confidence': top_tissue[1]['confidence'],
        'alternatives': [t[0] for t in alternatives if t[1]['composite_score'] > 0],
        'discriminators_present': discriminators_present,
        'discriminators_needed_to_confirm': discriminators_needed
    }
    
    return report
```

**Example output:**

```
CELL TYPE IDENTIFICATION RESULT:

Transcriptional Program: FOXI1+ Ion-Transporting Cell
Program Confidence: VERY HIGH (all methods concordant)

Most Likely Tissue: Pulmonary Ionocyte (Lung/Airway Epithelium)
Tissue Confidence: MEDIUM-HIGH

Evidence FOR lung ionocyte:
  ✅ All core genes present (FOXI1, V-ATPase, BSND)
  ✅ Literature confirms FOXI1 ionocytes in lung
  ✅ Fewer missing expected markers (1/4 missing)

Alternative Tissue: Kidney Intercalated Cell (Collecting Duct)
Tissue Confidence: MEDIUM

Evidence AGAINST kidney IC:
  ❌ Missing key discriminators: SLC4A1/AE1, SLC26A4, CAII
  ❌ These are expected in kidney ICs but absent

DISCRIMINATING GENES TO CONFIRM TISSUE:
  - CFTR (very high expression) → Confirms lung ionocyte
  - SLC4A1/AE1 OR SLC26A4 → Would indicate kidney IC
  - ATP6V0A4 + cochlear context → Would indicate inner ear

RECOMMENDATION:
  Search for CFTR in dataset. If CFTR is highly expressed, 
  confirms pulmonary ionocyte. If absent/low, reconsider kidney IC.
```

### Updated Comprehensive Workflow

The integrated workflow must now include tissue disambiguation:

```python
def identify_cell_type_comprehensive_v2(gene_list):
    """
    Comprehensive cell type identification with tissue disambiguation.
    """
    results = {
        'methods': {},
        'consensus': None,
        'confidence': None,
        'tissue_disambiguation': None
    }
    
    # Method 1: Check against known markers (fast)
    marker_result = check_known_markers(gene_list)
    results['methods']['markers'] = marker_result
    
    # NEW: Check for multi-tissue master TFs
    multitissue_tfs = check_multitissue_tf(gene_list)
    
    if multitissue_tfs:
        # TISSUE DISAMBIGUATION REQUIRED
        for tf_info in multitissue_tfs:
            tf_name = tf_info['tf']
            
            # Multi-tissue literature search
            lit_searches = multitissue_literature_search(tf_name, gene_list)
            lit_results = perform_literature_searches(lit_searches)
            
            # Check tissue discriminators
            disc_scores = check_tissue_discriminators(gene_list, tf_name)
            
            # Analyze missing markers
            missing = analyze_missing_markers(gene_list, tf_name)
            
            # Build confidence matrix
            tissue_confidence = build_tissue_confidence_matrix(
                gene_list, tf_name, lit_results, disc_scores, missing
            )
            
            results['tissue_disambiguation'] = {
                'tf': tf_name,
                'tissues': tf_info['tissues'],
                'confidence_matrix': tissue_confidence,
                'discriminators': disc_scores,
                'missing_analysis': missing
            }
    
    # Method 2: Enrichment analysis (high confidence)
    enrichment_result = enrichment_based_inference(gene_list)
    results['methods']['enrichment'] = enrichment_result
    
    # Method 3: GO/Pathway functional annotation
    functional_result = functional_annotation_inference(gene_list)
    results['methods']['functional'] = functional_result
    
    # Method 4: Expression database queries (for unknown genes)
    unknown_genes = identify_unknown_genes(gene_list, marker_result)
    if unknown_genes:
        expression_result = query_expression_databases_bulk(unknown_genes)
        results['methods']['expression_db'] = expression_result
    
    # Method 5: Co-expression analysis (if needed)
    if len(unknown_genes) > 0:
        coexpr_result = coexpression_inference(unknown_genes)
        results['methods']['coexpression'] = coexpr_result
    
    # Consensus building (now aware of tissue disambiguation)
    results['consensus'] = build_consensus_with_tissue_awareness(
        results['methods'],
        results.get('tissue_disambiguation')
    )
    
    results['confidence'] = assess_confidence_with_tissue_awareness(
        results['methods'],
        results.get('tissue_disambiguation')
    )
    
    return results
```

### Key Principles

1. **Recognize transcriptional programs accurately** (framework already does this well!)
2. **DO NOT assume single tissue** when master TF has multi-tissue expression
3. **Search ALL tissues** where master TF operates
4. **Use discriminating genes** to distinguish tissues
5. **Analyze missing markers** - what SHOULD be there but isn't?
6. **Report confidence appropriately:**
   - HIGH for transcriptional program
   - MEDIUM for tissue (if discriminators missing)
   - List alternatives and discriminators needed

### Success Metrics

**Before fix:**
- ❌ Reports: "Kidney IC, VERY HIGH confidence" (wrong tissue)

**After fix:**
- ✅ Reports: "FOXI1+ ion transporter, VERY HIGH confidence for program"
- ✅ Reports: "Most likely: Lung ionocyte (MEDIUM-HIGH)"
- ✅ Reports: "Alternative: Kidney IC (MEDIUM, but missing SLC4A1/SLC26A4)"
- ✅ Reports: "Discriminators needed: CFTR (lung) vs SLC4A1 (kidney)"

---

## Practical Example: Unknown Gene

**Scenario**: Gene "LRRC15" appears in your list but is not in the skill's marker database.

**Multi-Method Analysis:**

**Method 1: Known Markers**
- Result: NOT FOUND in predefined list

**Method 2: Enrichment Analysis**
```python
enrichr(['LRRC15'], 'Human_Gene_Atlas')
# Top result: "Fibroblasts" (p = 0.001)
```

**Method 3: GO/Pathway**
```python
go_enrichment(['LRRC15'])
# Top terms: "collagen fibril organization", "extracellular matrix"
# Inference: Likely fibroblast/stromal cell
```

**Method 4: Literature Mining**
```
PubMed search: "LRRC15 cell type"
# Abstracts mention: "cancer-associated fibroblast marker"
# Abstracts mention: "myofibroblast"
```

**Method 5: Expression Database**
```
Human Protein Atlas: LRRC15
# Cell type distribution: High in fibroblasts, CAFs
# Tissue: Connective tissue, stroma
```

**Method 6: Co-expression**
```
STRING: LRRC15 co-expressed with:
- COL1A1, FN1, POSTN (all fibroblast markers)
```

**Consensus:**
- **Cell Type**: Cancer-Associated Fibroblasts (CAFs)
- **Confidence**: HIGH
- **Evidence**: 5/5 methods concordant

---

## Special Cases and Challenges

### Case 1: Ambiguous Genes (Multiple Cell Types)

**Example**: CD44 (stem cells, some T cells, fibroblasts, tumor cells)

**Strategy:**
- Look at OTHER genes in the list
- Context matters: If list has tumor markers → likely tumor cells
- If list has T cell markers → likely T cells
- Check expression level patterns

### Case 2: Ubiquitous Genes

**Example**: GAPDH, ACTB (housekeeping genes)

**Strategy:**
- These don't indicate cell type
- Flag and exclude from analysis
- Focus on cell-type-restricted genes

### Case 3: Novel/Poorly Characterized Genes

**Example**: New lncRNAs, uncharacterized ORFs

**Strategy:**
- Rely heavily on co-expression
- Use guilt-by-association
- Lower confidence in final assignment
- May require experimental validation

### Case 4: Mixed Populations

**Example**: Gene list contains markers for multiple cell types

**Strategy:**
- Identify ALL represented cell types
- Report as "mixed population" or "TME signature"
- Estimate relative proportions if possible
- Consider biological context (tumor tissue = expected mixed)

---

## Recommended Default Workflow for Agents

```python
def agent_cell_type_identification(gene_list):
    """
    Recommended workflow for AI agent to identify cell types
    from ANY gene list, including unknown genes.
    """
    
    # STEP 1: Quick check against known markers (predefined in skill)
    print("[Step 1] Checking known marker databases...")
    known_result = check_marker_database(gene_list)
    
    if known_result['confidence'] == 'HIGH':
        # 3+ markers found, high confidence
        return known_result
    
    # STEP 2: Enrichment analysis (works for any genes)
    print("[Step 2] Running cell-type enrichment analysis...")
    enrichment_result = run_enrichment_analysis(
        gene_list,
        libraries=['Human_Gene_Atlas', 'ARCHS4_Tissues', 'Tabula_Sapiens']
    )
    
    # STEP 3: For genes not clearly identified, query additional resources
    unidentified_genes = identify_unknown_genes(gene_list, known_result)
    
    if len(unidentified_genes) > 0:
        print(f"[Step 3] {len(unidentified_genes)} genes need additional analysis...")
        
        # Functional annotation
        functional_result = go_pathway_analysis(unidentified_genes)
        
        # Web search for cell type associations (if needed)
        if len(unidentified_genes) <= 5:  # Limit to avoid too many queries
            web_results = {}
            for gene in unidentified_genes:
                web_results[gene] = web_search_cell_type(gene)
    
    # STEP 4: Build consensus
    print("[Step 4] Building consensus...")
    consensus = integrate_all_evidence(
        known_markers=known_result,
        enrichment=enrichment_result,
        functional=functional_result,
        web_evidence=web_results if 'web_results' in locals() else None
    )
    
    # STEP 5: Assess confidence and report
    final_result = {
        'cell_types': consensus['cell_types'],
        'confidence': consensus['confidence'],
        'evidence': {
            'known_markers': known_result,
            'enrichment': enrichment_result,
            'functional': functional_result if 'functional_result' in locals() else None
        },
        'method_concordance': consensus['concordance']
    }
    
    return final_result


def web_search_cell_type(gene):
    """
    Use web search to find cell type associations for a gene.
    """
    # Search pattern
    query = f"{gene} cell type expression marker"
    
    # This would use web_search tool
    results = web_search(query)
    
    # Extract cell type mentions from search results
    cell_type_mentions = extract_cell_types_from_text(results)
    
    return cell_type_mentions
```

---

## Key Recommendations for Agents

**Priority Order:**

1. **Always try enrichment analysis first** (works for any genes, statistical)
2. **Use GO/pathway annotation** as secondary evidence
3. **Web search for unknowns** (use sparingly, 3-5 genes max per query)
4. **Build consensus** from multiple methods (never rely on single method)
5. **Report confidence** explicitly (HIGH/MEDIUM/LOW)
6. **Show evidence** for transparency

**Never:**
- Claim you "don't recognize" genes without trying enrichment
- Rely only on predefined marker lists
- Give high confidence with single-method evidence
- Ignore genes just because they're unknown

**Always:**
- Try computational methods before giving up
- Combine multiple sources of evidence
- Report uncertainty honestly
- Provide reasoning for cell type assignment

---

## Summary

### The Agent Should:

✅ **First**: Check predefined markers (fast, high confidence if found)

✅ **Second**: Run enrichment analysis (works for ANY genes)

✅ **Third**: Use GO/pathway functional annotation

✅ **Fourth**: Query expression databases for unknown genes

✅ **Fifth**: Use web search strategically for critical unknowns

✅ **Always**: Build consensus from multiple methods

✅ **Always**: Report confidence level (HIGH/MEDIUM/LOW)

### The Agent Should NOT:

❌ Claim inability to identify cell types without trying multiple methods

❌ Rely solely on predefined marker lists

❌ Ignore genes not in skill database

❌ Give definitive answers without consensus

❌ Skip enrichment analysis (most powerful general method)

### Example Agent Response Pattern:

```
"I'll identify the cell types using multiple methods:

[Step 1] Checking against known marker databases...
- Found: CD8A → CD8+ T cells
- Not found: GENE_X, GENE_Y

[Step 2] Running enrichment analysis on all genes...
- Top enrichment: "Macrophages" (p = 1e-5)
- Supporting evidence: GENE_X, GENE_Y enriched in macrophage signatures

[Step 3] Functional annotation of unknown genes...
- GENE_X: GO term 'phagocytosis' → supports macrophage
- GENE_Y: Pathway 'TLR signaling' → supports macrophage

[Consensus]
Cell Type: Macrophages (likely M2-polarized based on specific markers)
Confidence: HIGH
Methods concordant: 3/3 (markers, enrichment, functional)
"
```

This approach ensures robust cell type identification regardless of whether genes are in predefined lists.

---

## Key References and Literature Support

### TME Phenotype Classification - Literature Support

The "Immune Paradox" / T-cell-excluded TLS-positive phenotype is well-documented in published literature:

#### Seminal Papers:

1. **Joshi NS, et al. (2015).** Regulatory T Cells in Tumor-Associated Tertiary Lymphoid Structures Suppress Anti-tumor T Cell Responses. *Immunity* 43(3):579-590.
   - **Key finding:** Tregs in TLS suppress anti-tumor responses in lung adenocarcinoma mouse model
   - **Mechanism:** Tregs suppress DC costimulatory ligand expression and T cell proliferation
   - **Result:** Treg depletion leads to tumor destruction
   - **Relevance:** Direct demonstration of TLS+ but immunosuppressed phenotype

2. **Tertiary Lymphoid Structure and CD8 T Cell Exclusion in Minimally Invasive Adenocarcinoma.** medRxiv. 2020.
   - **Key finding:** Early-stage lung adenocarcinoma with TLS present but CD8+ T cells excluded
   - **Composition:** Higher plasma cells, higher CD4+ Tregs, LOWER CD8+ cytotoxic T cells
   - **Mechanism:** Enrichment of follicular regulatory T cells in TLS follicles explains CD8+ exclusion
   - **Relevance:** Exact match to "Immune Paradox" phenotype

3. **Sautes-Fridman C, et al. (2019).** Tertiary lymphoid structures in the era of cancer immunotherapy. *Nature Reviews Cancer* 19(6):307-325.
   - **Key finding:** TLS heterogeneity - not all TLS are beneficial
   - **Composition:** TLS function depends on cellular composition (Tregs vs effector cells)
   - **Relevance:** Framework for understanding TLS variability

4. **Tertiary lymphoid structures in cancer: spatiotemporal heterogeneity, immune orchestration, and translational opportunities.** *Journal of Hematology & Oncology* (2025).
   - **Key finding:** TLS maturity stages - E-/PFL-TLS (immature) vs SFL-TLS (mature)
   - **Immature TLS:** Accompanied by Tregs, immunosuppressive, poor ICI response
   - **Mature TLS:** Functional germinal centers, anti-tumor
   - **Relevance:** Explains why TLS+ doesn't always predict good prognosis

5. **The pivotal role of tertiary lymphoid structures in the tumor immune microenvironment.** *Frontiers in Oncology* (2025).
   - **Key finding:** TLS have DUAL roles - anti-tumor vs pro-tumor depending on composition
   - **Tregs in TLS:** Predict worse outcomes despite TLS presence
   - **Mechanism:** Tregs suppress T cell function within organized TLS structure
   - **Relevance:** Validates "Immune Paradox" as recognized phenotype

6. **Regulatory T Cells: Barriers of Immune Infiltration Into the Tumor Microenvironment.** *Frontiers in Immunology* (2021).
   - **Key finding:** Tregs prevent T cell infiltration through multiple mechanisms
   - **CCL22/CCR4:** Treg recruitment pathway (targetable with mogamulizumab)
   - **Clinical data:** Mogamulizumab + nivolumab increased CD8+ T cells, decreased Tregs
   - **Relevance:** Therapeutic strategy for Treg depletion

#### Therapeutic Strategies - Evidence Base:

**IMpower150 Trial (FDA-approved):**
- Regimen: Atezolizumab + bevacizumab + carboplatin/paclitaxel
- Indication: Metastatic non-squamous NSCLC
- Mechanism: Anti-PD-L1 + vascular normalization + chemotherapy
- Relevance: Evidence-based combination therapy for lung adenocarcinoma

**Treg Depletion Studies:**
- Mogamulizumab (anti-CCR4): Phase I trials in solid tumors
- Metronomic cyclophosphamide: Preferentially depletes Tregs
- Combination approaches: Treg depletion + checkpoint blockade

**TLS Induction/Maturation:**
- Radiation + Treg therapy: Induces meningeal TLS formation
- STING agonists: Activate DCs, promote TLS maturation
- IL-7 adjuvants: Enhance TLS formation and function

### General TME Analysis References:

**Enrichment Databases:**
- PanglaoDB: Single-cell marker database
- CellMarker 2.0: Curated cell type markers
- Tabula Sapiens/Muris: Single-cell atlases
- Human Protein Atlas: Cell-type-specific expression

**Pathway Analysis:**
- Gene Ontology: Functional annotation
- KEGG: Pathway enrichment
- Reactome: Biological pathways
- MSigDB: Molecular signatures database

**TME Characterization:**
- TIMER2.0: Tumor immune estimation
- xCell: Cell type enrichment from expression data
- CIBERSORT: Immune cell deconvolution
- ImmuCellAI: Immune cell abundance identifier

### Key Insight from Literature:

**TLS are heterogeneous structures:**
- Mature TLS (SFL-TLS): Germinal centers, functional, anti-tumor
- Immature TLS (PFL-TLS): No germinal centers, may be immunosuppressive
- Treg-infiltrated TLS: Suppress anti-tumor responses despite structural organization
- Composition matters more than presence: Treg:CD8 ratio determines function

**Clinical Implications:**
- Don't assume TLS+ = good prognosis
- Assess TLS maturity (germinal centers?)
- Assess TLS composition (Tregs vs CD8+ ratio?)
- Assess vascular competence (can T cells enter?)
- Combination therapy addresses multiple barriers

---

## Version History

**Version 2.2 (Current):**
- Added comprehensive literature support for "Immune Paradox" phenotype
- Updated TME classification with evidence-based treatment recommendations
- Added references section with key papers and therapeutic strategies

**Version 2.1:**
- Added TME-aware DEG gene selection strategies
- Added stratified sampling approach for comprehensive TME analysis
- Added TME phenotype classification framework
- Added validation checklist and pitfall avoidance

**Version 2.0:**
- Added multi-tissue master TF disambiguation
- Added 6-method cell type identification approach
- Enhanced consensus building framework

**Version 1.0:**
- Initial framework for cell type identification beyond marker lists

