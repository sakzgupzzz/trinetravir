# HIV-1 biology — forward flag for Session 6B calibration interpretation

**Why this note exists (Schema v6 Part B4)**: the v1 training corpus (Wilk + Lee + Arunachalam + Schulte) contains exclusively acute RNA respiratory viruses (SARS-CoV-2 + Lee within-study IAV). The Lee 2025 HIV held-out cohort is the v1.5+ first retrovirus context. HIV biology differs from acute RNA respiratory viruses in three structural ways that affect the conserved-component transfer-learning question (Hypothesis H1):

## 1. Reverse transcription + integration (molecular)

HIV-1 is a retrovirus: viral genomic RNA is reverse-transcribed into double-stranded DNA, which integrates into the host genome via the viral integrase enzyme. This is fundamentally different from acute RNA viruses (SARS-CoV-2, IAV, RSV) which replicate cytoplasmically without genomic integration. The cellular sensors triggered by HIV are distinct: cGAS-STING (DNA sensing) is activated alongside MDA5/RIG-I (RNA sensing); the integrated provirus is a persistent host-genome modification that drives chronic transcriptional changes beyond the acute IFN response.

## 2. CD4 T cell tropism (cellular)

HIV-1 preferentially infects CD4 T cells (via CD4 receptor + CCR5/CXCR4 coreceptors). Monocytes can be infected at lower efficiency. CD8 T cells are NOT directly infected but expand clonally in response to HIV antigen. SARS-CoV-2 / IAV / RSV target respiratory epithelium primarily; PBMC response in those cohorts is the *bystander* response of circulating immune cells to systemic cytokine signals, not direct infection of PBMC.

For Session 6B Issue 30: the CD4T bucket is the primary biological focus for HIV. CD4T response in HIV is direct-infection biology (viral protein expression, integration-driven transcription). CD4T response in v1's SARS/IAV cohorts is bystander activation. These are different biologies — expected cross-context Pearson r on CD4T MVS subset is *very low* (predicted 0.00-0.20).

## 3. Chronic-by-definition (temporal)

Even "early" HIV (<6 months in Lee 2025) is chronic-by-definition: the integrated provirus is permanent. There is no acute-resolution phase. v1's SARS/IAV/RSV cohorts capture *acute* infection where the immune response is the immediate IFN-driven antiviral program. HIV early infection captures the *establishment of chronic persistent activation* — partially overlapping with acute IFN but dominated by chronic-state signatures (immune exhaustion markers, persistent activation, expansion of antigen-experienced T cells).

## Forward implications

Per METHODS_CHOICES Issue 30:
- **Expected outcome**: cross-context Pearson r on CD4T MVS subset is in [0.00, 0.20]. This is the *biologically defensible* outcome.
- **Surprising outcome**: cross-context r > 0.40 on CD4T MVS subset. Would require investigation — suggests v1's conserved component is capturing chronic-activation biology, not acute-IFN biology.
- **Concerning outcome**: cross-context r < -0.10 on CD4T MVS subset. Anti-correlation suggests HIV CD4T response is *opposite* to acute RNA virus CD4T response, which is also biologically interpretable (acute CD4T = transient activation + IFN; chronic HIV CD4T = exhaustion + immune dysregulation).

## References for paper write-up (Session 6B)

- Walker BD et al. HIV pathogenesis (review). *Nature* 580, 188-200 (2020).
- Yin K et al. CD4 T cells in HIV. *Annu Rev Immunol* 39, 525-547 (2021).
- Buggert M et al. T cell responses in early HIV. *J Exp Med* 217, e20200084 (2020).
- Doyle EH et al. HIV-1 and the type I IFN response. *Cell Host Microbe* 26, 5-9 (2019) — IFN-tone dysregulation in chronic HIV.
