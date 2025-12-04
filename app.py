import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from collections import Counter
import pandas as pd
import altair as alt

# ============================================
# CONFIGURACIÓN DE PÁGINA
# ============================================
st.set_page_config(
    page_title="Analizador de Secuencias ADN",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================
# CSS PERSONALIZADO PROFESIONAL
# ============================================
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono&display=swap');

    :root {
        --primary: #00D4AA;
        --secondary: #7C3AED;
        --accent: #F59E0B;
        --background: #0F172A;
        --surface: #1E293B;
        --text: #F1F5F9;
        --text-muted: #94A3B8;
    }

    .stApp {
        background: linear-gradient(135deg, #0F172A 0%, #1E293B 100%);
    }

    h1, h2, h3 {
        font-family: 'Inter', sans-serif !important;
        font-weight: 700 !important;
    }

    .main-header {
        background: linear-gradient(90deg, #00D4AA 0%, #7C3AED 50%, #F59E0B 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        font-size: 3rem;
        font-weight: 800;
        text-align: center;
        margin-bottom: 0.5rem;
        font-family: 'Inter', sans-serif;
    }

    .sub-header {
        text-align: center;
        color: #94A3B8;
        font-size: 1.1rem;
        margin-bottom: 2rem;
        font-family: 'Inter', sans-serif;
    }

    .metric-card {
        background: linear-gradient(145deg, #1E293B 0%, #334155 100%);
        border-radius: 16px;
        padding: 1.5rem;
        border: 1px solid rgba(255,255,255,0.1);
        box-shadow: 0 4px 20px rgba(0,0,0,0.3);
        transition: transform 0.3s ease, box-shadow 0.3s ease;
    }

    .metric-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 8px 30px rgba(0,212,170,0.2);
    }

    .metric-value {
        font-size: 2.5rem;
        font-weight: 700;
        background: linear-gradient(90deg, #00D4AA, #7C3AED);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-family: 'Inter', sans-serif;
    }

    .metric-label {
        color: #94A3B8;
        font-size: 0.9rem;
        text-transform: uppercase;
        letter-spacing: 1px;
        margin-bottom: 0.5rem;
        font-family: 'Inter', sans-serif;
    }

    .section-header {
        display: flex;
        align-items: center;
        gap: 0.75rem;
        margin: 2rem 0 1rem 0;
        padding-bottom: 0.5rem;
        border-bottom: 2px solid rgba(0,212,170,0.3);
    }

    .section-title {
        font-size: 1.5rem;
        font-weight: 600;
        color: #F1F5F9;
        font-family: 'Inter', sans-serif;
    }

    .sequence-box {
        background: #0F172A;
        border: 1px solid #334155;
        border-radius: 12px;
        padding: 1rem;
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.85rem;
        color: #00D4AA;
        overflow-x: auto;
        word-break: break-all;
        line-height: 1.6;
    }

    .sequence-box-protein {
        background: #0F172A;
        border: 1px solid #7C3AED;
        border-radius: 12px;
        padding: 1rem;
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.85rem;
        color: #F59E0B;
        overflow-x: auto;
        word-break: break-all;
    }

    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1E293B 0%, #0F172A 100%);
    }

    .stButton > button {
        background: linear-gradient(90deg, #00D4AA 0%, #7C3AED 100%);
        color: white;
        border: none;
        border-radius: 12px;
        padding: 0.75rem 2rem;
        font-weight: 600;
        font-size: 1rem;
        transition: all 0.3s ease;
        font-family: 'Inter', sans-serif;
    }

    .stButton > button:hover {
        transform: scale(1.02);
        box-shadow: 0 4px 20px rgba(0,212,170,0.4);
    }

    .stTextArea textarea {
        background: #0F172A;
        border: 2px solid #334155;
        border-radius: 12px;
        color: #00D4AA;
        font-family: 'JetBrains Mono', monospace;
    }

    .stTextArea textarea:focus {
        border-color: #00D4AA;
        box-shadow: 0 0 0 3px rgba(0,212,170,0.2);
    }

    .info-box {
        background: rgba(0,212,170,0.1);
        border-left: 4px solid #00D4AA;
        border-radius: 0 12px 12px 0;
        padding: 1rem;
        margin: 1rem 0;
    }

    .nuc-A { color: #FF6B6B; font-weight: bold; }
    .nuc-T { color: #4ECDC4; font-weight: bold; }
    .nuc-G { color: #45B7D1; font-weight: bold; }
    .nuc-C { color: #96CEB4; font-weight: bold; }
    .aa-stop { color: #EF4444; font-weight: bold; }

    .footer {
        text-align: center;
        padding: 2rem;
        color: #64748B;
        border-top: 1px solid #334155;
        margin-top: 3rem;
    }

    .footer a {
        color: #00D4AA;
        text-decoration: none;
    }

    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ============================================
# HEADER PRINCIPAL
# ============================================
st.markdown('<p class="main-header">Analizador de Secuencias ADN</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Herramienta  de análisis de secuencias de ADN para bioinformática</p>', unsafe_allow_html=True)

# ============================================
# SIDEBAR
# ============================================
with st.sidebar:
    st.markdown("## 🔬 Panel de Control")
    st.markdown("---")

    st.markdown("### 📊 Funcionalidades")
    st.markdown("""
    - ✅ Composición nucleotídica
    - ✅ Contenido GC/AT
    - ✅ Peso molecular
    - ✅ Transcripción (ADN → ARN)
    - ✅ Traducción (ARN → Proteína)
    - ✅ Detección de ORFs
    - ✅ Exportar resultados
    """)

    st.markdown("---")
    st.markdown("### ⚙️ Configuración")
    min_orf_length = st.slider("Longitud mínima ORF (bp)", 30, 150, 30, 15)

    st.markdown("---")
    st.markdown("### 📚 Referencias")
    st.markdown("""
    - [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
    - [Biopython Docs](https://biopython.org/)
    - [UniProt](https://www.uniprot.org/)
    """)

    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #64748B; font-size: 0.8rem;'>
        <p>Desarrollado para</p>
        <p><strong>Bioinformática - UNISON</strong></p>
        <p>2025</p>
    </div>
    """, unsafe_allow_html=True)

# ============================================
# INPUT DE SECUENCIA
# ============================================
st.markdown("""
<div class="section-header">
    <span style="font-size: 1.5rem;">📥</span>
    <span class="section-title">Entrada de Secuencia</span>
</div>
""", unsafe_allow_html=True)

col1, col2 = st.columns([3, 1])

with col1:
    ejemplo_secuencia = """ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
ATGCGTAAACTGAAACTTCTGCACGAATAATGCACCGACTAG
ATGGAAGCGAATCAATGCGCTGTACTCCTAG"""

    sequence_input = st.text_area(
        "Ingresa tu secuencia de ADN:",
        value=ejemplo_secuencia,
        height=150,
        label_visibility="collapsed",
        placeholder="Pega tu secuencia de ADN aquí..."
    )

with col2:
    st.markdown("<br>", unsafe_allow_html=True)
    st.info("💡 **Tip:** Puedes pegar secuencias en formato FASTA o texto plano.")

# Funciones auxiliares
def limpiar_secuencia(seq):
    lines = seq.strip().split("\n")
    cleaned = ""
    for line in lines:
        if not line.startswith(">"):
            cleaned += line.strip()
    return ''.join(c for c in cleaned.upper() if c in "ATGC")

def validar_secuencia(seq):
    return len(seq) > 0 and all(n in "ATGC" for n in seq)

def colorear_secuencia(seq):
    colored = ""
    for i, n in enumerate(seq):
        colored += f'<span class="nuc-{n}">{n}</span>'
        if (i + 1) % 10 == 0:
            colored += " "
        if (i + 1) % 60 == 0:
            colored += "<br>"
    return colored

def colorear_proteina(seq):
    colored = ""
    for i, aa in enumerate(seq):
        if aa == "*":
            colored += f'<span class="aa-stop">{aa}</span>'
        else:
            colored += aa
        if (i + 1) % 10 == 0:
            colored += " "
    return colored

sequence_clean = limpiar_secuencia(sequence_input)

# Botón de análisis
col1, col2, col3 = st.columns([1, 1, 1])
with col2:
    analyze_button = st.button("🔬 Analizar Secuencia", use_container_width=True, type="primary")

# ============================================
# ANÁLISIS
# ============================================
if analyze_button:
    if not validar_secuencia(sequence_clean):
        st.error("⚠️ Por favor ingresa una secuencia de ADN válida (solo A, T, G, C)")
    else:
        dna_seq = Seq(sequence_clean)

        nucleotide_counts = Counter(str(dna_seq))
        gc_content = gc_fraction(dna_seq) * 100
        at_content = 100 - gc_content

        try:
            peso_mol = molecular_weight(dna_seq, seq_type="DNA")
        except:
            peso_mol = len(dna_seq) * 330

        # MÉTRICAS BÁSICAS
        st.markdown("""
        <div class="section-header">
            <span style="font-size: 1.5rem;">📊</span>
            <span class="section-title">Métricas Básicas</span>
        </div>
        """, unsafe_allow_html=True)

        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">📏 Longitud</div>
                <div class="metric-value">{len(dna_seq):,}</div>
                <div style="color: #64748B; font-size: 0.9rem;">pares de bases</div>
            </div>
            """, unsafe_allow_html=True)

        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">🧪 Contenido GC</div>
                <div class="metric-value">{gc_content:.1f}%</div>
                <div style="color: #64748B; font-size: 0.9rem;">guanina + citosina</div>
            </div>
            """, unsafe_allow_html=True)

        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">🔬 Contenido AT</div>
                <div class="metric-value">{at_content:.1f}%</div>
                <div style="color: #64748B; font-size: 0.9rem;">adenina + timina</div>
            </div>
            """, unsafe_allow_html=True)

        with col4:
            st.markdown(f"""
            <div class="metric-card">
                <div class="metric-label">⚖️ Peso Molecular</div>
                <div class="metric-value">{peso_mol/1000:.1f}</div>
                <div style="color: #64748B; font-size: 0.9rem;">kDa (kilodaltons)</div>
            </div>
            """, unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        # COMPOSICIÓN DE NUCLEÓTIDOS
        st.markdown("""
        <div class="section-header">
            <span style="font-size: 1.5rem;">🧬</span>
            <span class="section-title">Composición de Nucleótidos</span>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns([1, 2])

        with col1:
            names = {"A": "Adenina", "T": "Timina", "G": "Guanina", "C": "Citosina"}
            nuc_data = []
            for nuc in "ATGC":
                count = nucleotide_counts.get(nuc, 0)
                pct = count / len(dna_seq) * 100
                nuc_data.append({
                    "Nucleótido": f"{names[nuc]} ({nuc})",
                    "Conteo": count,
                    "Porcentaje": f"{pct:.1f}%"
                })

            df_nuc = pd.DataFrame(nuc_data)
            st.dataframe(df_nuc, hide_index=True, use_container_width=True)

            st.markdown("""
            <div class="info-box">
                <strong>📐 Reglas de Chargaff:</strong><br>
                A ≈ T y G ≈ C en ADN de doble cadena
            </div>
            """, unsafe_allow_html=True)

            ratio_at = nucleotide_counts.get('A', 0) / max(nucleotide_counts.get('T', 1), 1)
            ratio_gc = nucleotide_counts.get('G', 0) / max(nucleotide_counts.get('C', 1), 1)
            st.markdown(f"**Ratio A/T:** {ratio_at:.2f} | **Ratio G/C:** {ratio_gc:.2f}")

        with col2:
            df_chart = pd.DataFrame({
                "Nucleótido": list("ATGC"),
                "Conteo": [nucleotide_counts.get(n, 0) for n in "ATGC"],
                "Nombre": ["Adenina", "Timina", "Guanina", "Citosina"]
            })

            bars = alt.Chart(df_chart).mark_bar(
                cornerRadiusTopLeft=8,
                cornerRadiusTopRight=8
            ).encode(
                x=alt.X("Nucleótido:N", axis=alt.Axis(labelAngle=0, labelFontSize=14)),
                y=alt.Y("Conteo:Q", axis=alt.Axis(labelFontSize=12)),
                color=alt.Color("Nucleótido:N", 
                    scale=alt.Scale(domain=["A", "T", "G", "C"], 
                                   range=["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"]),
                    legend=None),
                tooltip=["Nombre:N", "Conteo:Q"]
            ).properties(height=300)

            text = bars.mark_text(
                align='center', baseline='bottom', dy=-5, fontSize=14, fontWeight='bold'
            ).encode(text='Conteo:Q', color=alt.value('#F1F5F9'))

            st.altair_chart((bars + text).configure_axis(grid=False).configure_view(strokeWidth=0), use_container_width=True)

            pie = alt.Chart(df_chart).mark_arc(innerRadius=50).encode(
                theta=alt.Theta("Conteo:Q"),
                color=alt.Color("Nucleótido:N",
                    scale=alt.Scale(domain=["A", "T", "G", "C"],
                                   range=["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"])),
                tooltip=["Nombre:N", "Conteo:Q"]
            ).properties(height=200)
            st.altair_chart(pie, use_container_width=True)

        # TRANSCRIPCIÓN Y TRADUCCIÓN
        st.markdown("""
        <div class="section-header">
            <span style="font-size: 1.5rem;">🔄</span>
            <span class="section-title">Transcripción y Traducción</span>
        </div>
        """, unsafe_allow_html=True)

        mrna = dna_seq.transcribe()
        protein = dna_seq.translate()
        complement = dna_seq.complement()
        reverse_comp = dna_seq.reverse_complement()

        tab1, tab2, tab3 = st.tabs(["🧬 Secuencias de ADN", "📋 ARN Mensajero", "🧫 Proteína"])

        with tab1:
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**ADN Original (5' → 3')**")
                st.markdown(f'<div class="sequence-box">{colorear_secuencia(str(dna_seq))}</div>', unsafe_allow_html=True)
            with col2:
                st.markdown("**Cadena Complementaria (3' → 5')**")
                st.markdown(f'<div class="sequence-box">{colorear_secuencia(str(complement))}</div>', unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("**Reverso Complementario (5' → 3')**")
            st.markdown(f'<div class="sequence-box">{colorear_secuencia(str(reverse_comp))}</div>', unsafe_allow_html=True)

        with tab2:
            st.markdown("**ARN mensajero (mRNA)**")
            mrna_str = str(mrna)
            mrna_colored = mrna_str.replace("A", '<span style="color:#FF6B6B">A</span>').replace("U", '<span style="color:#4ECDC4">U</span>').replace("G", '<span style="color:#45B7D1">G</span>').replace("C", '<span style="color:#96CEB4">C</span>')
            st.markdown(f'<div class="sequence-box">{mrna_colored}</div>', unsafe_allow_html=True)

            st.markdown("""
            <div class="info-box">
                <strong>📝 Proceso:</strong> Durante la transcripción, la Timina (T) se reemplaza por Uracilo (U)
            </div>
            """, unsafe_allow_html=True)

        with tab3:
            st.markdown("**Secuencia de Aminoácidos**")
            protein_str = str(protein)
            st.markdown(f'<div class="sequence-box-protein">{colorear_proteina(protein_str)}</div>', unsafe_allow_html=True)

            col1, col2, col3 = st.columns(3)
            col1.metric("Aminoácidos totales", len(protein))
            col2.metric("Codones STOP", protein_str.count("*"))
            col3.metric("Met (inicio)", protein_str.count("M"))

        # ORFs
        st.markdown("""
        <div class="section-header">
            <span style="font-size: 1.5rem;">🔍</span>
            <span class="section-title">Marcos de Lectura Abiertos (ORFs)</span>
        </div>
        """, unsafe_allow_html=True)

        def encontrar_orfs(sequence, min_length=30):
            orfs = []
            seq_str = str(sequence)
            for frame in range(3):
                i = frame
                while i < len(seq_str) - 2:
                    if seq_str[i:i+3] == "ATG":
                        for j in range(i + 3, len(seq_str) - 2, 3):
                            stop = seq_str[j:j+3]
                            if stop in ["TAA", "TAG", "TGA"]:
                                orf_len = j + 3 - i
                                if orf_len >= min_length:
                                    orf_seq = seq_str[i:j+3]
                                    orf_protein = str(Seq(orf_seq).translate())
                                    orfs.append({
                                        "Marco": f"+{frame + 1}",
                                        "Inicio": i + 1,
                                        "Fin": j + 3,
                                        "Longitud (bp)": orf_len,
                                        "Aminoácidos": len(orf_protein),
                                        "Proteína": orf_protein if len(orf_protein) <= 20 else orf_protein[:20] + "..."
                                    })
                                i = j + 3
                                break
                        else:
                            i += 3
                    else:
                        i += 3
            return orfs

        orfs = encontrar_orfs(dna_seq, min_orf_length)

        if orfs:
            df_orfs = pd.DataFrame(orfs)
            st.dataframe(df_orfs, hide_index=True, use_container_width=True)

            orf_chart = alt.Chart(df_orfs).mark_bar(height=20).encode(
                x=alt.X("Inicio:Q", scale=alt.Scale(domain=[0, len(dna_seq)])),
                x2="Fin:Q",
                y=alt.Y("Marco:N"),
                color=alt.Color("Marco:N", scale=alt.Scale(scheme="category10")),
                tooltip=["Marco", "Inicio", "Fin", "Longitud (bp)"]
            ).properties(height=150, title="Mapa de ORFs")
            st.altair_chart(orf_chart, use_container_width=True)
        else:
            st.info(f"No se encontraron ORFs con longitud mínima de {min_orf_length} bp")

        # EXPORTAR
        st.markdown("""
        <div class="section-header">
            <span style="font-size: 1.5rem;">💾</span>
            <span class="section-title">Exportar Resultados</span>
        </div>
        """, unsafe_allow_html=True)

        col1, col2, col3 = st.columns(3)

        reporte = f"""DNA SEQUENCE ANALYZER - REPORTE
================================
Longitud: {len(dna_seq)} bp
Contenido GC: {gc_content:.2f}%
Contenido AT: {at_content:.2f}%
Peso Molecular: {peso_mol:.0f} Da

COMPOSICIÓN:
A: {nucleotide_counts.get('A', 0)} ({nucleotide_counts.get('A', 0)/len(dna_seq)*100:.1f}%)
T: {nucleotide_counts.get('T', 0)} ({nucleotide_counts.get('T', 0)/len(dna_seq)*100:.1f}%)
G: {nucleotide_counts.get('G', 0)} ({nucleotide_counts.get('G', 0)/len(dna_seq)*100:.1f}%)
C: {nucleotide_counts.get('C', 0)} ({nucleotide_counts.get('C', 0)/len(dna_seq)*100:.1f}%)

SECUENCIAS:
ADN: {str(dna_seq)}
ARN: {str(mrna)}
Proteína: {str(protein)}

ORFs encontrados: {len(orfs)}
"""

        with col1:
            st.download_button("📄 Descargar TXT", reporte, "reporte.txt", "text/plain", use_container_width=True)

        with col2:
            fasta = f">Secuencia length={len(dna_seq)}\n{str(dna_seq)}\n>mRNA\n{str(mrna)}\n>Protein\n{str(protein)}"
            st.download_button("🧬 Descargar FASTA", fasta, "sequences.fasta", "text/plain", use_container_width=True)

        with col3:
            st.download_button("📊 Descargar CSV", df_nuc.to_csv(index=False), "composicion.csv", "text/csv", use_container_width=True)

# FOOTER
st.markdown("""
<div class="footer">
    <p><strong>Analizador de Secuencias ADN</strong></p>
    <p>Desarrollado con Python, Streamlit y Biopython</p>
    <p>Clase de Bioinformática - Universidad de Sonora | Maestro: Aaron Lara Ordoñez</p>
</div>
""", unsafe_allow_html=True)