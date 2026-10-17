# 映射规则定义文件####
# 文件: Code/mapping_rules.R
# 说明: 定义癌种三层分类映射规则和共用筛选规则
# 使用: 在主分析脚本中通过 source("Code/mapping_rules.R") 调用
#
# 内容包括：
#   1. 恶性/良性肿瘤筛选正则表达式（三个数据库共用）
#   2. 三层分类体系映射规则
#
# 三层分类体系:
#   - 第一层 (Primary Organ/Tissue): 原发器官/组织
#     - 具体器官（如Lung, Kidney, Breast等）
#     - Rare（稀有癌种，如肾上腺皮质癌）
#     - Unknown（无法确定原发灶）
#
#   - 第二层 (Histologic Category): 在第一层基础上内部划分
#     - 例如: Lung → SCLC, NSCLC, Unknown
#     - 例如: Skin → Melanoma, SCC, BCC, Unknown
#
#   - 第三层 (TCGA Entity): 在第二层基础上细分到TCGA代码
#     - 例如: NSCLC → LUAD, LUSC, Unknown
#     - 例如: Melanoma → SKCM, Unknown

message("加载映射规则...")

# 共用筛选规则定义（FAERS/JADER两个数据库共用）####
## 恶性肿瘤筛选规则####
### A. 强恶性词（组织学/疾病名，高特异）
primary_strong <- "(?i)\\b(ADENOCARCINOMA|SARCOMA|MELANOMA|LYMPHOMA|LEUKA?EMI[AE]|MYELOMA|GLIOBLASTOMA|\\w+BLASTOMA|MESOTHELIOMA|CHOLANGIOCARCINOMA|HODGKIN|MYELODYSPLASTIC|MYELOPROLIFERATIVE|KAPOSI|MYCOSIS\\s+FUNGOIDES|PLASMACYTOMA)\\b"

### B. 泛化组织学恶性（不指明原发，但仍强提示恶性）
generalized_histology <- "(?i)\\b(CARCINOMA|MALIGNANT\\s+(TUMOU?R|NEOPLASM|DISEASE))\\b"
### C. 泛化恶性诊断（未特指部位/癌种）
generalized_malignant <- "(?i)\\b(CANCER|MALIGNANCY|NEOPLASM\\s+MALIGNANT|MALIGNANT\\s+(NEOPLASM|TUMOU?R))\\b"

### D. 转移（作为恶性确认的一部分；同时也单独用于 metastasis_only 判定）
metastasis_pattern <- "(?i)\\b(METASTAS\\w*|METASTATIC|SECONDARY\\s+(MALIGNANT|CANCER|CARCINOMA|NEOPLASM))\\b"

### E. 转移到某部位（用于排除器官原发癌误判，如"METASTATIC IN LIVER"不应映射为LIHC）
met_to_site_pattern <- "(?i)\\b(METASTA[ST]\\w*|SECONDARY)\\b.{0,30}\\b(IN|TO)\\b"

### 合并：确定恶性肿瘤总pattern
malignant_confirmed_pattern <- paste0(
  "(", primary_strong, "|", generalized_histology, "|", generalized_malignant, "|", metastasis_pattern, ")"
)

## 明确良性肿瘤筛选规则####
# 包括：BENIGN标记、良性肿瘤后缀、以及特定的良性/交界性肿瘤名称
# 2026-01-18更新：新增HAEMANGIOBLASTOMA(WHO Grade I)、AMELOBLASTOMA、CHONDROBLASTOMA
benign_pattern <- "(?i)\\b(BENIGN|NON[- ]?MALIGNANT|NONMALIGNANT|ADENOMA|LIPOMA|FIBROMA|HA?EMANGIOMA|PAPILLOMA|HAMARTOMA|NEVUS|NAEVUS|LEIOMYOMA|MYOMA|LIPOBLASTOMA|MYELOBLASTOMA|ANGIOMYOFIBROBLASTOMA|HA?EMANGIOBLASTOMA|AMELOBLASTOMA|CHONDROBLASTOMA)\\b"

## 假阳性排除规则（被误识别为恶性的非肿瘤诊断）####
# 这些PT包含恶性肿瘤关键词但实际上不是癌症诊断
# 例如：Metastatic Pulmonary Embolism（转移性肺栓塞）是并发症，不是癌症
# 例如：Anaplastic Lymphoma Kinase Gene Mutation（ALK基因突变）是基因特征，不是癌症诊断
false_positive_pattern <- "(?i)\\b(METASTATIC\\s+PULMONARY\\s+EMBOLISM|PULMONARY\\s+EMBOLISM\\s+METASTATIC|ANAPLASTIC\\s+LYMPHOMA\\s+KINASE\\s+GENE\\s+MUTATION|ANAPLASTIC\\s+LYMPHOMA\\s+KINASE\\s+GENE\\s+AND\\s+NUCLEOPHOSMIN\\s+GENE\\s+FUSION\\s+OVEREXPRESSION)\\b"

## 良恶性不明肿瘤筛选规则####
# 注：THYMOMA/PHEOCHROMOCYTOMA/PARAGANGLIOMA已移至癌种映射（THYM/PCPG）
uncertain_pattern <- "(?i)\\b(NEOPLASM|TUMOU?R|MENINGIOMA|TERATOMA)\\b"

## 泛化恶性（原发不可判定）筛选规则####
# 这些PT确认恶性但无法判定原发部位，不应进入主分析Others
# 使用宽松匹配，覆盖NOS/UNSPECIFIED/UNKNOWN PRIMARY等变体
# 注意：SARCOMA已移除，因为可以映射到TCGA的SARC代码
generalized_only_pattern <- "(?i)^(NEOPLASM MALIGNANT|MALIGNANT NEOPLASM|CANCER|MALIGNANCY|CARCINOMA|ADENOCARCINOMA|SQUAMOUS CELL CARCINOMA|LYMPHOMA|CANCER PAIN|MALIGNANT DISEASE|TUMOU?R MALIGNANT|CARCINOMA OF UNKNOWN PRIMARY|UNKNOWN PRIMARY|MALIGNANT NEOPLASM NOS|CANCER NOS|CARCINOMA NOS|NEOPLASM NOS|MALIGNANCY NOS|CANCER UNSPECIFIED|MALIGNANT TUMOU?R|DISSEMINATED CANCER|ADVANCED CANCER|TERMINAL CANCER|STAGE IV CANCER|METASTATIC CANCER|METASTATIC DISEASE)(\\s*(NOS|UNSPECIFIED|,|\\s))*$"

# 癌种映射函数（PT → cancer_type，TCGA缩写）####
# 说明：将适应证PT映射到TCGA癌种代码
# 返回：TCGA代码（如ACC、BLCA等）或NA
map_cancer_type <- function(pt_text) {
  pt <- toupper(pt_text)

  # ACC - 肾上腺皮质癌
  if (grepl("ADRENOCORTICAL CARCINOMA", pt)) {
    return("ACC")
  }

  # BLCA - 膀胱尿路上皮癌
  if (
    grepl("\\b(BLADDER|BLCA)\\b", pt) ||
      (grepl("\\bUROTHELIAL\\b", pt) && !grepl("URETER|RENAL\\s+PELVIS|UPPER\\s+TRACT|UTUC|URETHR", pt))
  ) {
    if (grepl("CANCER|CARCINOMA|NEOPLASM|TUMOU?R", pt) &&
      !grepl("BENIGN|GALLBLADDER", pt)) {
      return("BLCA")
    }
  }

  # BRCA - 乳腺癌
  if (grepl("BREAST|BRCA", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR|ADENOCARCINOMA|ADENOSQUAMOUS|SQUAMOUS", pt) &&
    !grepl("BENIGN", pt)) {
    return("BRCA")
  }

  # CESC - 宫颈癌
  if (grepl("CERVICAL|CERVIX", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN", pt)) {
    return("CESC")
  }

  # CHOL - 胆管癌
  if (grepl("CHOLANGIOCARCINOMA|BILE DUCT|\\bBILIARY\\b|HEPATOBILIARY|PANCREATICOBILIARY", pt) &&
    grepl("CANCER|CARCINOMA|ADENOCARCINOMA|NEOPLASM|TUMOU?R|MALIGNANT", pt) &&
    !grepl("BENIGN|PANCREAS|PANCREATIC", pt)) {
    return("CHOL")
  }


  # COAD - 结肠腺癌
  if (grepl("\\b(COLON|COAD)\\b", pt) &&
    !grepl("\\b(COLORECTAL|RECTOSIGMOID)\\b", pt) &&
    grepl("\\b(ADENOCARCINOMA|ADENOSQUAMOUS|CANCER|CARCINOMA|NEOPLASM|TUMOU?R|MALIGNANT)\\b", pt) &&
    !grepl("\\bBENIGN\\b", pt)) {
    return("COAD")
  }

  # DLBC - 弥漫性大B细胞淋巴瘤
  if (grepl("DIFFUSE LARGE B-CELL", pt) &&
    grepl("LYMPHOMA", pt) &&
    !grepl("BENIGN", pt)) {
    return("DLBC")
  }

  # ESCA - 食管癌
  if (grepl("\\b(ESOPHAGEAL|OESOPHAGEAL)\\b", pt) &&
    !grepl("\\b(GASTRO[- ]?(ESOPHAGEAL|OESOPHAGEAL)|GEJ|GASTROESOPHAGEAL|GASTROOESOPHAGEAL)\\b", pt) &&
    grepl("\\b(CANCER|CARCINOMA|NEOPLASM|TUMOU?R|ADENOCARCINOMA|SQUAMOUS|ADENOSQUAMOUS)\\b", pt) &&
    !grepl("\\bBENIGN\\b", pt)) {
    return("ESCA")
  }

  # GBM - 胶质母细胞瘤
  if (grepl("GLIOBLASTOMA|GBM", pt) && !grepl("BENIGN", pt)) {
    return("GBM")
  }

  # HNSC - 头颈部鳞状细胞癌
  if (grepl("HEAD AND NECK|HEAD|NECK|HNSC|LIP|ORAL|PALATE|OROPHARYNGEAL|TONGUE|HYPOPHARYNGEAL|NASOPHARYNX|NASOPHARYNGEAL|PHARYNGEAL|PHARYNX|GLOTTIS|LARYNGEAL|NASAL CAVITY|NASAL|TONSIL", pt) &&
    grepl("SQUAMOUS|CANCER|CARCINOMA|NEOPLASM|TUMOU?R", pt) &&
    !grepl("BENIGN|LUNG|ESOPHAGEAL|OESOPHAGEAL|SALIVARY|PAROTID|SUBMANDIBULAR|SUBLINGUAL|ADENOID CYSTIC|MUCOEPIDERMOID|SALIVARY DUCT|ACINIC", pt)) {
    return("HNSC")
  }

  # KICH - 嫌色细胞癌
  if (grepl("CHROMOPHOBE|KICH", pt) && grepl("RENAL|KIDNEY", pt) && !grepl("BENIGN", pt)) {
    return("KICH")
  }

  # KIRC - 肾透明细胞癌
  if (grepl("CLEAR CELL", pt) && grepl("KIDNEY|RENAL", pt) && !grepl("BENIGN|NON", pt)) {
    return("KIRC")
  }

  # KIRP - 肾乳头状细胞癌
  if (grepl("PAPILLARY|KIRP", pt) && grepl("RENAL|KIDNEY", pt) && !grepl("BENIGN", pt)) {
    return("KIRP")
  }

  # LAML - 急性髓系白血病
  if (grepl("ACUTE|LAML", pt) &&
    grepl("ERYTHROID|MEGAKARYOCYTIC|MONOCYTIC|MYELOID|MYELOMONOCYTIC|PROMYELOCYTIC", pt) &&
    grepl("LEUKA?EMIA", pt) && !grepl("BENIGN", pt)) {
    return("LAML")
  }

  # LGG - 低级别胶质瘤（需要LOW GRADE关键词）
  if (grepl("ASTROCYTOMA|LGG|GLIOMA|GLIOBLASTOMA|GLIOMATOSIS|EPENDYMOMA|OLIGODENDROGLIOMA|MEDULLOBLASTOMA|OLIGOASTROCYTOMA", pt) &&
    grepl("LOW GRADE|LOW[- ]GRADE", pt) &&
    !grepl("BENIGN", pt)) {
    return("LGG")
  }

  # LIHC - 肝细胞癌
  if (grepl("HEPATOCELLULAR|LIVER|HEPATIC|LIHC", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN|CHOLANGIO|BILE\\s+DUCT|BILIARY|HEPATOBILIARY|PANCREATICOBILIARY|MIXED\\s+HEPATOCELLULAR", pt) &&
    !grepl(met_to_site_pattern, pt, perl = TRUE)) {
    return("LIHC")
  }

  # LUAD - 肺腺癌
  if (grepl("LUNG|LUAD", pt) && grepl("ADENOCARCINOMA|ADENOSQUAMOUS", pt) &&
    !grepl("BENIGN", pt) && !grepl(met_to_site_pattern, pt, perl = TRUE)) {
    return("LUAD")
  }

  # LUSC - 肺鳞状细胞癌
  if (grepl("LUNG|LUSC", pt) && grepl("SQUAMOUS", pt) &&
    !grepl("BENIGN", pt) && !grepl(met_to_site_pattern, pt, perl = TRUE)) {
    return("LUSC")
  }

  # 注：无法细分为LUAD/LUSC的肺癌归入Others（不使用通用LUNG）

  # MESO - 间皮瘤
  if (grepl("MESOTHELIOMA|MESO", pt) &&
    !grepl("BENIGN|ACROMESOMELIC DYSPLASIA|MESOBLASTIC NEPHROMA|MESOTHERAPY", pt)) {
    return("MESO")
  }

  # OV - 卵巢癌
  if (grepl("\\b(OVARIAN|OVARY|OV)\\b", pt) &&
    grepl("\\b(CANCER|CARCINOMA|NEOPLASM|TUMOU?R)\\b", pt) &&
    !grepl("\\bBENIGN\\b", pt) &&
    !grepl("\\b(GERM\\s*CELL|EMBRYONAL|CHORIOCARCINOMA|YOLK\\s*SAC|DYSGERMINOMA|TERATOMA)\\b", pt) &&
    !grepl("\\b(STROMAL|SEX\\s*CORD|GRANULOSA|SERTOLI|LEYDIG|THECOMA|FIBROMA)\\b", pt) &&
    !grepl("\\b(REMOVAL|RESECTION|EXCISION|SURGERY|BIOPSY)\\b", pt)
  ) {
    return("OV")
  }

  # PAAD - 胰腺癌
  if (grepl("PANCREATIC|PANCREAS|PAAD", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN", pt) && !grepl(met_to_site_pattern, pt, perl = TRUE)) {
    return("PAAD")
  }

  # PCPG - 嗜铬细胞瘤和副神经节瘤
  if (grepl("PHEOCHROMOCYTOMA|PHAEOCHROMOCYTOMA|PARAGANGLIOMA|PARAGANGLION(\\s+NEOPLASM)?|PCPG|GLOMUS JUGULARE TUMOU?R", pt) &&
    !grepl("BENIGN|PSEUDOPHAEOCHROMOCYTOMA|GLOMANGIOPERICYTOMA|\\bGLOMUS\\s+TUMOU?R\\b", pt)) {
    return("PCPG")
  }

  # PRAD - 前列腺癌
  if (grepl("PROSTATE|PRAD", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN", pt)) {
    return("PRAD")
  }

  # READ - 直肠腺癌
  if (grepl("\\b(RECTAL|RECTUM|READ)\\b", pt) &&
    !grepl("\\b(COLORECTAL|RECTOSIGMOID)\\b", pt) &&
    grepl("\\b(ADENOCARCINOMA|ADENOSQUAMOUS|CANCER|CARCINOMA|NEOPLASM|TUMOU?R|MALIGNANT)\\b", pt) &&
    !grepl("\\bBENIGN\\b", pt)) {
    return("READ")
  }

  # SARC - 软组织肉瘤
  if (grepl("SOFT TISSUE SARCOMA|LIPOSARCOMA|LEIOMYOSARCOMA|DEDIFFERENTIATED LIPOSARCOMA|UNDIFFERENTIATED PLEOMORPHIC SARCOMA|MYXOFIBROSARCOMA|MALIGNANT PERIPHERAL NERVE SHEATH TUMORS|MPNST|SARCOMA", pt) &&
    !grepl("BENIGN|UTERINE|KAPOSI|OSTEOSARCOMA|CHONDROSARCOMA|GIST|GASTROINTESTINAL STROMAL|RHABDOMYOSARCOMA|EWING|UTERUS|SARCOMATOID\\s+CARCINOMA|EXCISION|RESECTION|SURGERY|PSEUDOSARCOMA", pt)) {
    return("SARC")
  }
  # UVM - 葡萄膜黑色素瘤
  if (grepl("UVEAL|UVM|OCULAR|INTRAOCULAR|CHOROID|CILIARY\\s*BODY|IRIS|CONJUNCTIV", pt) &&
    grepl("MELANOMA", pt) &&
    !grepl("BENIGN", pt) &&
    !grepl("ANTI-?MELANOMA|ANTIBODY|DERMATOMYOSITIS|MDA5|MDA-5", pt)) {
    return("UVM")
  }

  # SKCM - 皮肤黑色素瘤
  if (grepl("MELANOMA|LENTIGO\\s+MALIGNA|MALIGNANT\\s+BLUE\\s+NAEVUS", pt) &&
    !grepl("UVEAL|UVM|OCULAR|INTRAOCULAR|CHOROID|CILIARY\\s*BODY|IRIS|CONJUNCTIV", pt) &&
    !grepl("GASTROINTESTINAL|GENITOURINARY|URETHRAL|OVARIAN|PULMONARY|CENTRAL\\s+NERVOUS\\s+SYSTEM|CNS", pt) &&
    !grepl("ANTI-?MELANOMA|ANTIBODY|DERMATOMYOSITIS|MDA5|MDA-5", pt) &&
    !grepl("BENIGN", pt)) {
    return("SKCM")
  }

  # STAD - 胃腺癌
  if (grepl("STOMACH|GASTRIC|STAD", pt) &&
    grepl("ADENOCARCINOMA|CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN|\\bADENOMA\\b|OVARY|OVARIAN|CYSTADENOCARCINOMA\\s+OVARY|PANCREATIC\\s+CYSTADENOMA|LYMPHOMATOSUM|PANCREATIC", pt)) {
    return("STAD")
  }

  # TGCT - 睾丸生殖细胞肿瘤
  if (grepl("TESTICULAR|TESTIS", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN", pt)) {
    return("TGCT")
  }

  # THCA - 甲状腺癌
  if (grepl("THYROID|THCA", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR", pt) &&
    !grepl("BENIGN", pt)) {
    return("THCA")
  }

  # THYM - 胸腺瘤（包括THYMOMA和THYMIC CARCINOMA）
  if (grepl("THYMOMA|THYMIC CARCINOMA|THYMIC CANCER|NEOPLASM OF THYMUS", pt) &&
    !grepl("BENIGN", pt)) {
    return("THYM")
  }
  # UCS - 子宫癌肉瘤
  if (grepl("CARCINOSARCOMA|MIXED\\s+MULLERIAN|MMMT|M\\u00dcLLER|MUELLER|LEIOMYOSARCOMA|SARCOMA\\s+UTERUS|UTERINE\\s+SARCOMA", pt) &&
    !grepl("BENIGN", pt)) {
    return("UCS")
  }
  # UCEC - 子宫内膜癌
  if (grepl("UTERINE|UTERUS|ENDOMETRIAL|UCEC", pt) &&
    grepl("CANCER|CARCINOMA|NEOPLASM|TUMOR|ADENOCARCINOMA|ADENOSQUAMOUS|SQUAMOUS", pt) &&
    !grepl("BENIGN|ADENOMA|IN\\s*SITU|STAGE\\s*0|ADNEXA|CERVIX|CERVICAL|SARCOMA|CARCINOSARCOMA|MIXED\\s+MULLERIAN|MMMT|OVAR|FALLOPIAN|TUBAL|SALPING", pt)
  ) {
    return("UCEC")
  }
  return(NA_character_)
}

# TCGA assigned_type 到器官的映射（优先级最高）####
# 格式：tcga_to_organ <- list("TCGA_CODE" = "Organ Name", ...)
# 说明：当 assigned_type 有值且非 "Others" 时，优先使用此映射
# 基于论文图b的分类体系

tcga_to_organ <- list(
  # Head/neck（头颈部）
  "HNSC" = "Head/neck",   # 头颈部鳞状细胞癌
  
  # Pancreas（胰腺）
  "PAAD" = "Pancreas",    # 胰腺癌
  
  # Biliary tract（胆道）
  "CHOL" = "Biliary tract", # 胆管癌
  
  # Liver（肝）
  "LIHC" = "Liver",       # 肝细胞癌
  
  # Esophagus/stomach（食管/胃）
  "STAD" = "Esophagus/stomach", # 胃腺癌
  "ESCA" = "Esophagus/stomach", # 食管癌
  
  # Bowel（肠道）- 包括结直肠、小肠、阑尾等
  "COAD" = "Bowel",       # 结肠腺癌
  "READ" = "Bowel",       # 直肠腺癌
  
  # Peritoneum（腹膜）
  # 无对应TCGA代码
  
  # Ovary（卵巢）
  "OV" = "Ovary",         # 卵巢癌
  
  # Breast（乳腺）
  "BRCA" = "Breast",      # 乳腺癌
  
  # Uterus（子宫）
  "UCEC" = "Uterus",      # 子宫内膜癌
  "UCS" = "Uterus",       # 子宫癌肉瘤
  
  # Cervix（宫颈）
  "CESC" = "Cervix",      # 宫颈癌
  
  # Prostate（前列腺）
  "PRAD" = "Prostate",    # 前列腺癌
  
  # Kidney（肾）
  "KIRC" = "Kidney",      # 肾透明细胞癌
  "KIRP" = "Kidney",      # 肾乳头状细胞癌
  "KICH" = "Kidney",      # 嫌色细胞癌
  
  # Bladder（膀胱）
  "BLCA" = "Bladder",     # 膀胱尿路上皮癌
  
  # Lung（肺）
  "LUAD" = "Lung",        # 肺腺癌
  "LUSC" = "Lung",        # 肺鳞状细胞癌
  
  # Pleura（胸膜）
  "MESO" = "Pleura",      # 间皮瘤
  
  # Thyroid（甲状腺）
  "THCA" = "Thyroid",     # 甲状腺癌
  
  # Thymus（胸腺）
  "THYM" = "Thymus",      # 胸腺瘤
  
  # Soft tissue（软组织）- 包括GIST、软组织肉瘤
  "SARC" = "Soft tissue", # 软组织肉瘤
  
  # Bone（骨）
  # 无对应TCGA代码
  
  # Brain（脑）
  "GBM" = "Brain",        # 胶质母细胞瘤
  "LGG" = "Brain",        # 低级别胶质瘤
  
  # PNS（外周神经系统）
  "PCPG" = "PNS",         # 嗜铬细胞瘤和副神经节瘤
  
  # Skin（皮肤）
  "SKCM" = "Skin",        # 皮肤黑色素瘤
  
  # Hematologic（血液系统）- 非实体瘤
  "LAML" = "Hematologic", # 急性髓系白血病
  "DLBC" = "Hematologic", # 弥漫性大B细胞淋巴瘤
  
  # Eye（眼）
  "UVM" = "Eye",          # 葡萄膜黑色素瘤
  
  # Other/Rare（其他/罕见）- CUP等
  "ACC" = "Rare",        # 肾上腺皮质癌
  "TGCT" = "Rare"        # 睾丸生殖细胞肿瘤
)

# 第二层映射规则：Histologic Category（器官内部分类）####
# 说明：在第一层器官的基础上进行内部分类
# 每个器官有自己的第二层分类体系
# 格式：organ_histology_xxx <- list(分类名 = c(PT1, PT2, ...))

## 第二层分类定义（每个器官的内部分类）####

### Lung（肺）第二层分类####
# SCLC: 小细胞肺癌
# NSCLC: 非小细胞肺癌
# Unknown: 只标注肺癌，无法确定类型
lung_histology <- list(
  "SCLC" = c(
    "Small Cell Lung Cancer",
    "Small Cell Lung Cancer Extensive",
    "Small Cell Lung Cancer Limited",
    "Small Cell Lung Cancer Recurrent",
    "Small Cell Lung Cancer Metastatic",
    "Small Cell Lung Cancer Extensive Stage"
  ),
  "NSCLC" = c(
    "Non-Small Cell Lung Cancer",
    "Non-Small Cell Lung Cancer Recurrent",
    "Non-Small Cell Lung Cancer Metastatic",
    "Her2 Mutant Non-Small Cell Lung Cancer",
    # 腺癌
    "Lung Adenocarcinoma Metastatic",
    "Bronchioloalveolar Carcinoma",
    # 鳞癌
    "Lung Squamous Cell Carcinoma Metastatic",
    # 大细胞癌
    "Large Cell Lung Cancer",
    "Large Cell Lung Cancer Recurrent",
    "Large Cell Lung Cancer Metastatic",
    # 肉瘤样癌
    "Sarcomatoid Carcinoma Of The Lung"
  ),
  "Unknown/Others" = c(
    "Lung Neoplasm Malignant",
    "Lung Carcinoma Cell Type Unspecified",
    "Lung Carcinoma Cell Type Unspecified Recurrent",
    "Bronchial Carcinoma",
    "Lung Cancer Metastatic",
    "Metastatic Bronchial Carcinoma",
    "Respiratory Tract Carcinoma In Situ",
    "Neuroendocrine Tumour Of The Lung Metastatic",
    "Primary Pulmonary Melanoma"
  )
)

### Kidney（肾）第二层分类####
# RCC: 肾细胞癌
# TCC: 尿路上皮癌（肾盂/输尿管）
# Unknown/Others: 只标注肾癌，无法确定类型
kidney_histology <- list(
  "RCC" = c(
    "Renal Cell Carcinoma",
    "Metastatic Renal Cell Carcinoma",
    "Renal Cell Carcinoma Recurrent",
    "Non-Renal Cell Carcinoma Of Kidney",
    "Hereditary Leiomyomatosis Renal Cell Carcinoma"
  ),
  "TCC" = c(
    "Transitional Cell Cancer Of The Renal Pelvis And Ureter",
    "Transitional Cell Cancer Of Renal Pelvis And Ureter Metastatic",
    "Transitional Cell Cancer Of The Renal Pelvis And Ureter Recurrent",
    "Transitional Cell Cancer Of The Renal Pelvis And Ureter Localised",
    "Transitional Cell Cancer Of The Renal Pelvis And Ureter Regional"
  ),
  "Unknown/Others" = c(
    "Renal Cancer",
    "Renal Cancer Recurrent",
    "Renal Cancer Metastatic",
    "Malignant Neoplasm Of Renal Pelvis",
    "Nephroblastoma"
  )
)

### Bladder（膀胱）第二层分类####
# Urothelial: 尿路上皮癌
# Unknown/Others: 只标注膀胱癌
bladder_histology <- list(
  "Urothelial" = c(
    "Transitional Cell Carcinoma",
    "Transitional Cell Carcinoma Urethra",
    "Transitional Cell Carcinoma Recurrent",
    "Transitional Cell Carcinoma Metastatic",
    "Bladder Transitional Cell Carcinoma Metastatic",
    "Urinary Tract Carcinoma In Situ"
  ),
  "Unknown/Others" = c(
    "Urethral Cancer",
    "Ureteric Cancer",
    "Ureteric Cancer Recurrent",
    "Ureteric Cancer Local",
    "Ureteric Cancer Metastatic",
    "Metastatic Carcinoma Of The Bladder",
    "Urethral Cancer Metastatic",
    "Malignant Urinary Tract Neoplasm Metastatic"
  )
)

### Skin（皮肤）第二层分类####
# Melanoma: 黑色素瘤
# SCC: 鳞状细胞癌
# BCC: 基底细胞癌
# Unknown/Others: 其他或无法确定
skin_histology <- list(
  "Melanoma" = c(
    "Metastatic Malignant Melanoma",
    "Urethral Melanoma Metastatic"
  ),
  "SCC" = c(
    "Squamous Cell Carcinoma Of Skin",
    "Skin Squamous Cell Carcinoma Recurrent",
    "Skin Squamous Cell Carcinoma Metastatic"
  ),
  "BCC" = c(
    "Basal Cell Carcinoma",
    "Basal Cell Carcinoma Metastatic",
    "Basosquamous Carcinoma",
    "Basosquamous Carcinoma Of Skin"
  ),
  "Unknown/Others" = c(
    "Skin Cancer",
    "Carcinoma In Situ Of Skin",
    "Skin Cancer Metastatic",
    "Neuroendocrine Carcinoma Of The Skin",
    "Cutaneous Lymphoma",
    "Cutaneous B-Cell Lymphoma",
    "Sebaceous Carcinoma",
    "Eccrine Carcinoma",
    "Pilomatrix Carcinoma",
    "Trichoblastic Carcinoma"
  )
)

### Brain（脑）第二层分类####
# Glioma: 胶质瘤（包括GBM和LGG）
# Other: 其他脑肿瘤
brain_histology <- list(
  "Glioma" = c(
    "Brain Neoplasm Malignant",
    "Astroblastoma",
    "Metastatic Glioma"
  ),
  "Medulloblastoma" = c(
    "Medulloblastoma",
    "Medulloblastoma Recurrent"
  ),
  "Other" = c(
    "Central Nervous System Lymphoma",
    "Aesthesioneuroblastoma",
    "Haemangioblastoma",
    "Central Nervous System Melanoma",
    "Central Nervous System Leukaemia",
    "Choroid Plexus Carcinoma",
    "Pinealoblastoma",
    "Pineal Parenchymal Neoplasm Malignant",
    "Central Nervous System Neuroblastoma",
    "Primitive Neuroectodermal Tumour Metastatic",
    "Malignant Meningioma Metastatic",
    # 原Unknown合并至此
    "Brain Cancer Metastatic",
    "Malignant Neoplasm Of Spinal Cord",
    "Pituitary Neoplasm Malignant Recurrent",
    "Pituitary Neoplasm Malignant",
    "Pituitary Cancer Metastatic"
  )
)

### Hematologic（血液系统）第二层分类####
# Lymphoma: 淋巴瘤
# Leukemia: 白血病
# Plasma_Cell: 浆细胞肿瘤
# MDS_MPN: 骨髓增生异常/骨髓增殖性肿瘤
# Unknown: 无法确定类型
hematologic_histology <- list(
  "Lymphoma" = c(
    # B细胞淋巴瘤
    "B-Cell Lymphoma", "B-Cell Lymphoma Refractory", "B-Cell Lymphoma Recurrent",
    "B-Cell Unclassifiable Lymphoma High Grade", "B-Cell Unclassifiable Lymphoma Low Grade",
    "High-Grade B-Cell Lymphoma", "High Grade B-Cell Lymphoma Refractory",
    "Primary Mediastinal Large B-Cell Lymphoma",
    "B-Cell Small Lymphocytic Lymphoma",
    # 套细胞淋巴瘤
    "Mantle Cell Lymphoma", "Mantle Cell Lymphoma Recurrent", "Mantle Cell Lymphoma Refractory",
    # 滤泡性淋巴瘤
    "Follicular Lymphoma", "Follicular Lymphoma Refractory", "Follicular Lymphoma Recurrent",
    # T细胞淋巴瘤
    "T-Cell Lymphoma", "T-Cell Lymphoma Refractory", "T-Cell Lymphoma Recurrent",
    "Cutaneous T-Cell Lymphoma", "Peripheral T-Cell Lymphoma Unspecified",
    "Angioimmunoblastic T-Cell Lymphoma",
    # 霍奇金淋巴瘤
    "Hodgkin's Disease", "Hodgkin's Disease Nodular Sclerosis",
    "Hodgkin's Disease Refractory", "Hodgkin's Disease Recurrent",
    # 非霍奇金淋巴瘤
    "Non-Hodgkin's Lymphoma", "Non-Hodgkin's Lymphoma Recurrent", "Non-Hodgkin's Lymphoma Refractory",
    # 间变性大细胞淋巴瘤
    "Anaplastic Large-Cell Lymphoma",
    # 边缘区淋巴瘤
    "Marginal Zone Lymphoma", "Marginal Zone Lymphoma Recurrent",
    # 伯基特淋巴瘤
    "Burkitt's Lymphoma", "Burkitt's Lymphoma Recurrent",
    # 其他淋巴瘤
    "Lymphoma", "Lymphoma Transformation", "Mycosis Fungoides"
  ),
  "Leukemia" = c(
    # 慢性淋巴细胞白血病
    "Chronic Lymphocytic Leukaemia", "Chronic Lymphocytic Leukaemia Recurrent",
    # 慢性髓系白血病
    "Chronic Myeloid Leukaemia", "Chronic Myeloid Leukaemia Recurrent",
    # 急性淋巴细胞白血病
    "Acute Lymphocytic Leukaemia", "Acute Lymphocytic Leukaemia Recurrent",
    # 其他白血病
    "Leukaemia", "Leukaemia Recurrent", "Acute Leukaemia", "Chronic Leukaemia",
    "Hairy Cell Leukaemia", "Natural Killer-Cell Leukaemia"
  ),
  "Plasma_Cell" = c(
    "Plasma Cell Myeloma", "Plasma Cell Myeloma In Remission",
    "Plasma Cell Myeloma Recurrent", "Plasma Cell Myeloma Refractory",
    "Plasmacytoma", "Plasma Cell Leukaemia"
  ),
  "MDS_MPN" = c(
    "Myelodysplastic Syndrome", "Myeloproliferative Neoplasm",
    "Myelodysplastic Syndrome With Excess Blasts",
    "Chronic Myelomonocytic Leukaemia"
  ),
  "Unknown/Others" = c(
    "Haematological Malignancy"
  )
)

### Head_neck（头颈部）第二层分类####
# SCC: 鳞状细胞癌
# Salivary: 唾液腺癌
# Unknown/Others: 无法确定类型
head_neck_histology <- list(
  "SCC" = c(
    "Squamous Cell Carcinoma Of The Parotid Gland"
  ),
  "Salivary" = c(
    "Salivary Gland Cancer", "Her2 Positive Salivary Gland Cancer",
    "Mucoepidermoid Carcinoma", "Adenoid Cystic Carcinoma",
    "Adenoid Cystic Carcinoma Of Salivary Gland",
    "Adenocarcinoma Of Salivary Gland",
    "Acinic Cell Carcinoma Of Salivary Gland",
    "Salivary Gland Cancer Recurrent", "Metastatic Salivary Gland Cancer",
    "Mucoepidermoid Carcinoma Of Salivary Gland",
    "Carcinoma Ex-Pleomorphic Adenoma",
    "Myoepithelial Carcinoma Of Salivary Gland"
  ),
  "Unknown/Others" = c(
    "Throat Cancer", "Gingival Cancer",
    "Epiglottic Cancer", "Epiglottic Carcinoma",
    "Tracheal Cancer", "Postcricoid Cancer",
    "Ear Neoplasm Malignant", "External Ear Neoplasm Malignant",
    "Malignant Neoplasm Of Auricular Cartilage",
    "Head And Neck Cancer Metastatic",
    "Tongue Cancer Metastatic", "Sinus Cancer Metastatic"
  )
)

### Esophagus_stomach（食管胃）第二层分类####
# SCC: 鳞状细胞癌（主要食管）
# Adenocarcinoma: 腺癌
# Unknown/Others: 无法确定类型
esophagus_stomach_histology <- list(
  "SCC" = c(
    "Oesophageal Squamous Cell Carcinoma Metastatic"
  ),
  "Adenocarcinoma" = c(
    # 无直接匹配PT，通过TCGA映射
  ),
  "Unknown/Others" = c(
    "Gastrooesophageal Cancer", "Gastrooesophageal Cancer Recurrent",
    "Oesophageal Cancer Metastatic", "Metastatic Gastric Cancer",
    "Neuroendocrine Carcinoma Of The Oesophagus"
  )
)

### Bowel（肠道）第二层分类####
# Adenocarcinoma: 腺癌
# Unknown/Others: 无法确定类型
bowel_histology <- list(
  "Adenocarcinoma" = c(
    "Colorectal Adenocarcinoma", "Intestinal Adenocarcinoma",
    "Small Intestine Adenocarcinoma", "Adenocarcinoma Of Appendix",
    "Gastrointestinal Adenocarcinoma", "Mucinous Adenocarcinoma Of Appendix"
  ),
  "SCC" = c(
    "Anal Squamous Cell Carcinoma"
  ),
  "Unknown/Others" = c(
    "Colorectal Cancer", "Rectosigmoid Cancer", "Anal Cancer", "Appendix Cancer",
    "Colorectal Cancer Metastatic", "Colorectal Cancer Recurrent",
    "Small Intestine Carcinoma", "Gastrointestinal Carcinoma",
    "Colorectal Carcinoma", "Colon Cancer Metastatic", "Rectal Cancer Metastatic"
  )
)

### Liver（肝）第二层分类####
# HCC: 肝细胞癌
# Cholangiocarcinoma: 胆管细胞癌
# Unknown/Others: 无法确定类型
liver_histology <- list(
  "HCC" = c(
    "Hepatic Cancer Metastatic"
  ),
  "Hepatoblastoma" = c(
    "Hepatoblastoma", "Hepatoblastoma Recurrent"
  ),
  "Unknown/Others" = c(
    # 大多数肝癌PT需要通过TCGA映射
  )
)

### Breast（乳腺）第二层分类####
# Carcinoma: 乳腺癌（TCGA只有BRCA一个分类）
# Lymphoma: 乳腺淋巴瘤
# Unknown/Others: 无法确定类型
breast_histology <- list(
  "Carcinoma" = c(
    "Breast Cancer Metastatic"
  ),
  "Lymphoma" = c(
    "Primary Breast Lymphoma",
    "Breast Implant-Associated Anaplastic Large Cell Lymphoma"
  ),
  "Unknown/Others" = c(
    # 其他需要通过TCGA映射
  )
)

### Ovary（卵巢）第二层分类####
# Epithelial: 上皮性卵巢癌
# Germ_Cell: 生殖细胞肿瘤
# Stromal: 性索间质肿瘤
# Unknown/Others: 无法确定类型
ovary_histology <- list(
  "Epithelial" = c(
    "Ovarian Cancer Metastatic"
  ),
  "Germ_Cell" = c(
    "Ovarian Germ Cell Cancer", "Ovarian Germ Cell Embryonal Carcinoma",
    "Ovarian Embryonal Carcinoma"
  ),
  "Stromal" = c(
    "Ovarian Stromal Cancer"
  ),
  "Unknown/Others" = c(
    "Fallopian Tube Cancer", "Fallopian Tube Cancer Metastatic"
  )
)

### Uterus（子宫）第二层分类####
# Endometrial: 子宫内膜癌
# Sarcoma: 子宫肉瘤
# Choriocarcinoma: 绒毛膜癌
# Unknown/Others: 无法确定类型
uterus_histology <- list(
  "Endometrial" = c(
    "Endometrial Cancer", "Endometrial Cancer Metastatic"
  ),
  "Sarcoma" = c(
    # UCS通过TCGA映射
  ),
  "Choriocarcinoma" = c(
    "Choriocarcinoma", "Choriocarcinoma Recurrent", "Metastatic Choriocarcinoma"
  ),
  "Unknown/Others" = c(
    "Malignant Neoplasm Of Uterine Adnexa",
    "Uterine Carcinoma In Situ", "Metastatic Uterine Cancer"
  )
)

### Soft_tissue（软组织）第二层分类####
# Sarcoma: 肉瘤
# GIST: 胃肠道间质瘤
# Unknown/Others: 无法确定类型
soft_tissue_histology <- list(
  "Sarcoma" = c(
    "Ewing's Sarcoma", "Ewing's Sarcoma Recurrent", "Ewing-Like Sarcoma",
    "Extra-Osseous Ewing's Sarcoma", "Ewing's Sarcoma Metastatic",
    "Kaposi's Sarcoma", "Kaposi's Sarcoma Aids Related",
    "Sarcoma Metastatic", "Liposarcoma Metastatic", "Leiomyosarcoma Metastatic",
    "Angiosarcoma Metastatic"
  ),
  "GIST" = c(
    "Gastrointestinal Stromal Cancer"
  ),
  "Unknown/Others" = c(
    "Soft Tissue Cancer", "Retroperitoneal Cancer", "Retroperitoneal Neoplasm Metastatic"
  )
)

# 第三层映射规则：TCGA Entity（TCGA代码）####
# 说明：在第二层的基础上细分到TCGA代码
# 格式：tcga_entity_xxx <- list(TCGA_CODE = c(PT1, PT2, ...), ...)

## TCGA assigned_type 到器官的映射（保留用于第一层映射）####
tcga_to_organ <- list(
  # Head/neck（头颈部）
  "HNSC" = "Head/neck",   # 头颈部鳞状细胞癌
  # Pancreas（胰腺）
  "PAAD" = "Pancreas",    # 胰腺癌
  # Biliary tract（胆道）
  "CHOL" = "Biliary tract", # 胆管癌
  # Liver（肝）
  "LIHC" = "Liver",       # 肝细胞癌
  # Esophagus/stomach（食管/胃）
  "STAD" = "Esophagus/stomach", # 胃腺癌
  "ESCA" = "Esophagus/stomach", # 食管癌
  # Bowel（肠道）
  "COAD" = "Bowel",       # 结肠腺癌
  "READ" = "Bowel",       # 直肠腺癌
  # Ovary（卵巢）
  "OV" = "Ovary",         # 卵巢癌
  # Breast（乳腺）
  "BRCA" = "Breast",      # 乳腺癌
  # Uterus（子宫）
  "UCEC" = "Uterus",      # 子宫内膜癌
  "UCS" = "Uterus",       # 子宫癌肉瘤
  # Cervix（宫颈）
  "CESC" = "Cervix",      # 宫颈癌
  # Prostate（前列腺）
  "PRAD" = "Prostate",    # 前列腺癌
  # Kidney（肾）
  "KIRC" = "Kidney",      # 肾透明细胞癌
  "KIRP" = "Kidney",      # 肾乳头状细胞癌
  "KICH" = "Kidney",      # 嫌色细胞癌
  # Bladder（膀胱）
  "BLCA" = "Bladder",     # 膀胱尿路上皮癌
  # Lung（肺）
  "LUAD" = "Lung",        # 肺腺癌
  "LUSC" = "Lung",        # 肺鳞状细胞癌
  # Pleura（胸膜）
  "MESO" = "Pleura",      # 间皮瘤
  # Thyroid（甲状腺）
  "THCA" = "Thyroid",     # 甲状腺癌
  # Thymus（胸腺）
  "THYM" = "Thymus",      # 胸腺瘤
  # Soft tissue（软组织）
  "SARC" = "Soft tissue", # 软组织肉瘤
  # Brain（脑）
  "GBM" = "Brain",        # 胶质母细胞瘤
  "LGG" = "Brain",        # 低级别胶质瘤
  # PNS（外周神经系统）
  "PCPG" = "PNS",         # 嗜铬细胞瘤和副神经节瘤
  # Skin（皮肤）
  "SKCM" = "Skin",        # 皮肤黑色素瘤
  # Hematologic（血液系统）
  "LAML" = "Hematologic", # 急性髓系白血病
  "DLBC" = "Hematologic", # 弥漫性大B细胞淋巴瘤
  # Eye（眼）
  "UVM" = "Eye",          # 葡萄膜黑色素瘤
  # Other/Rare（其他/罕见）
  "ACC" = "Rare",         # 肾上腺皮质癌
  "TGCT" = "Rare"         # 睾丸生殖细胞肿瘤（并入Rare）
)

## TCGA代码到第二层分类的映射####
# 2026-01-19 根据用户校准版本更新，与three_layer_mapping保持一致
tcga_to_histology <- list(
  # Lung第二层
  "LUAD" = "NSCLC",       # 肺腺癌 → 非小细胞肺癌
  "LUSC" = "NSCLC",       # 肺鳞癌 → 非小细胞肺癌
  
  # Kidney第二层
  "KIRC" = "Renal cell carcinoma",   # 肾透明细胞癌
  "KIRP" = "Renal cell carcinoma",   # 肾乳头状细胞癌
  "KICH" = "Renal cell carcinoma",   # 嫌色细胞癌
  
  # Bladder第二层
  "BLCA" = "Urothelial carcinoma",   # 膀胱尿路上皮癌
  
  # Skin第二层
  "SKCM" = "Melanoma",    # 皮肤黑色素瘤 → 黑色素瘤
  
  # Brain第二层
  "GBM" = "Glioma",       # 胶质母细胞瘤 → 胶质瘤
  "LGG" = "Glioma",       # 低级别胶质瘤 → 胶质瘤
  
  # Hematologic第二层
  "DLBC" = "Lymphoma",    # 弥漫性大B细胞淋巴瘤 → 淋巴瘤
  "LAML" = "Leukemia",    # 急性髓系白血病 → 白血病
  
  # Head/neck第二层
  "HNSC" = "Head/neck SCC",   # 头颈部鳞癌
  
  # Esophagus/stomach第二层
  "ESCA" = "Esophageal Cancer",   # 食管癌
  "STAD" = "Stomach Cancer",      # 胃癌
  
  # Bowel第二层
  "COAD" = "Colorectal cancer",   # 结肠癌 → 结直肠癌
  "READ" = "Colorectal cancer",   # 直肠癌 → 结直肠癌
  
  # Liver第二层
  "LIHC" = "HCC",         # 肝细胞癌
  
  # Breast第二层
  "BRCA" = "Breast Cancer",   # 乳腺癌
  
  # Ovary第二层
  "OV" = "Ovary cancer",      # 卵巢癌
  
  # Uterus第二层
  "UCEC" = "Endometrial",     # 子宫内膜癌
  "UCS" = "Sarcoma",          # 子宫癌肉瘤 → 肉瘤
  
  # Cervix第二层
  "CESC" = "Cervical cancer", # 宫颈癌
  
  # Prostate第二层
  "PRAD" = "Prostate cancer", # 前列腺癌
  
  # Pancreas第二层
  "PAAD" = "Pancreatic cancer", # 胰腺癌
  
  # Biliary tract第二层
  "CHOL" = "Cholangiocarcinoma", # 胆管癌
  
  # Thyroid第二层
  "THCA" = "Thyroid Cancer",  # 甲状腺癌
  
  # Thymus第二层
  "THYM" = "Thymoma",     # 胸腺瘤
  
  # Pleura第二层
  "MESO" = "Mesothelioma", # 间皮瘤
  
  # Soft tissue第二层
  "SARC" = "Soft tissue sarcoma", # 软组织肉瘤
  
  # PNS第二层
  "PCPG" = "PCPG",        # 嗜铬细胞瘤/副神经节瘤
  
  # Eye第二层
  "UVM" = "Uveal Melanoma", # 葡萄膜黑色素瘤
  
  # Testis第二层
  "TGCT" = "Testicular Cancer", # 睾丸癌
  
  # Rare第二层
  "ACC" = "ACC"           # 肾上腺皮质癌
)

## TCGA第三层（NSCLC内部细分）####
nsclc_tcga <- list(
  "LUAD" = c(
    "Lung Adenocarcinoma Metastatic",
    "Bronchioloalveolar Carcinoma"
  ),
  "LUSC" = c(
    "Lung Squamous Cell Carcinoma Metastatic"
  ),
  "Unknown/Others" = c(
    "Non-Small Cell Lung Cancer",
    "Non-Small Cell Lung Cancer Recurrent",
    "Non-Small Cell Lung Cancer Metastatic",
    "Her2 Mutant Non-Small Cell Lung Cancer",
    "Large Cell Lung Cancer",
    "Large Cell Lung Cancer Recurrent",
    "Large Cell Lung Cancer Metastatic",
    "Sarcomatoid Carcinoma Of The Lung"
  )
)

## TCGA第三层（RCC内部细分）####
rcc_tcga <- list(
  "KIRC" = c(
    # 透明细胞癌，多数通过TCGA assigned_type映射
  ),
  "KIRP" = c(
    # 乳头状细胞癌，多数通过TCGA assigned_type映射
  ),
  "KICH" = c(
    # 嫌色细胞癌，多数通过TCGA assigned_type映射
  ),
  "Unknown/Others" = c(
    "Renal Cell Carcinoma",
    "Metastatic Renal Cell Carcinoma",
    "Renal Cell Carcinoma Recurrent",
    "Non-Renal Cell Carcinoma Of Kidney",
    "Hereditary Leiomyomatosis Renal Cell Carcinoma"
  )
)

## TCGA第三层（Glioma内部细分）####
glioma_tcga <- list(
  "GBM" = c(
    # 胶质母细胞瘤，多数通过TCGA assigned_type映射
  ),
  "LGG" = c(
    # 低级别胶质瘤，多数通过TCGA assigned_type映射
  ),
  "Unknown/Others" = c(
    "Brain Neoplasm Malignant",
    "Astroblastoma",
    "Metastatic Glioma"
  )
)

## PT名称到组织学分类的精细映射规则####
## 说明：基于PT名称中的关键词进行组织学分类
## 优先级：PT关键词 > TCGA代码 > 器官默认分类

### 腺癌相关PT####
histology_Adenocarcinoma <- c(
  # 明确包含adenocarcinoma的PT
  "Colorectal Adenocarcinoma",
  "Intestinal Adenocarcinoma",
  "Gallbladder Adenocarcinoma",
  "Small Intestine Adenocarcinoma",
  "Gastrointestinal Adenocarcinoma",
  "Adenocarcinoma Of Appendix",
  "Mucinous Adenocarcinoma Of Appendix",
  "Adenocarcinoma Of Salivary Gland",
  "Adenocarcinoma",
  "Adenocarcinoma Metastatic",
  "Vaginal Adenocarcinoma",
  "Vulvar Adenocarcinoma",
  # 肺腺癌
  "Lung Adenocarcinoma Metastatic",
  "Bronchioloalveolar Carcinoma"  # 细支气管肺泡癌（现归类为原位腺癌）
)

### 鳞状细胞癌相关PT####
histology_Squamous_Cell_Carcinoma <- c(
  "Squamous Cell Carcinoma",
  "Squamous Cell Carcinoma Of Skin",
  "Squamous Cell Carcinoma Of The Vulva",
  "Squamous Cell Carcinoma Of The Vagina",
  "Squamous Cell Carcinoma Of The Parotid Gland",
  "Anal Squamous Cell Carcinoma",
  "Gallbladder Squamous Cell Carcinoma",
  "Metastatic Squamous Cell Carcinoma",
  "Skin Squamous Cell Carcinoma Recurrent",
  "Skin Squamous Cell Carcinoma Metastatic",
  "Penile Squamous Cell Carcinoma",
  # 肺鳞癌
  "Lung Squamous Cell Carcinoma Metastatic",
  # 食管鳞癌
  "Oesophageal Squamous Cell Carcinoma Metastatic",
  # 腺鳞癌（含鳞癌成分）
  "Adenosquamous Cell Carcinoma",
  "Basosquamous Carcinoma",
  "Basosquamous Carcinoma Of Skin"
)

### 尿路上皮癌相关PT####
histology_Urothelial_Carcinoma <- c(
  "Transitional Cell Carcinoma",
  "Transitional Cell Carcinoma Urethra",
  "Transitional Cell Carcinoma Recurrent",
  "Transitional Cell Carcinoma Metastatic",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter",
  "Transitional Cell Cancer Of Renal Pelvis And Ureter Metastatic",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Recurrent",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Localised",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Regional",
  "Bladder Transitional Cell Carcinoma Metastatic",
  "Urinary Tract Carcinoma In Situ"
)

### 小细胞癌相关PT####
histology_Small_Cell_Carcinoma <- c(
  "Small Cell Lung Cancer",
  "Small Cell Lung Cancer Extensive",
  "Small Cell Lung Cancer Limited",
  "Small Cell Lung Cancer Recurrent",
  "Small Cell Lung Cancer Metastatic",
  "Small Cell Lung Cancer Extensive Stage",
  "Small Cell Carcinoma"
)

### 黑色素瘤相关PT####
histology_Melanoma <- c(
  "Metastatic Malignant Melanoma",
  "Primary Pulmonary Melanoma",
  "Central Nervous System Melanoma",
  "Gastrointestinal Melanoma",
  "Genitourinary Melanoma",
  "Urethral Melanoma Metastatic"
)

### 肉瘤相关PT####
histology_Sarcoma <- c(
  # Ewing肉瘤
  "Ewing's Sarcoma",
  "Ewing's Sarcoma Recurrent",
  "Ewing-Like Sarcoma",
  "Extra-Osseous Ewing's Sarcoma",
  "Ewing's Sarcoma Metastatic",
  "Extra-Osseous Ewing's Sarcoma Metastatic",
  # Kaposi肉瘤
  "Kaposi's Sarcoma",
  "Kaposi's Sarcoma Aids Related",
  "Kaposi's Sarcoma Classical Type",
  "Kaposi Sarcoma Inflammatory Cytokine Syndrome",
  "Immune Reconstitution Inflammatory Syndrome Associated Kaposi's Sarcoma",
  # 其他肉瘤
  "Sarcomatoid Carcinoma",
  "Sarcomatoid Carcinoma Of The Lung",
  "Sarcoma Metastatic",
  "Liposarcoma Metastatic",
  "Leiomyosarcoma Metastatic",
  "Angiosarcoma Metastatic",
  "Malignant Fibrous Histiocytoma Metastatic",
  "Malignant Haemangiopericytoma Metastatic",
  "Malignant Mesenchymoma Metastatic",
  # 骨肉瘤
  "Osteosarcoma Metastatic",
  "Osteosarcoma Recurrent",
  "Extraskeletal Osteosarcoma Metastatic",
  # 软骨肉瘤
  "Chondrosarcoma Metastatic",
  "Extraskeletal Chondrosarcoma Metastatic",
  # GIST
  "Gastrointestinal Stromal Cancer"
)

### 胶质瘤相关PT####
histology_Glioma <- c(
  "Brain Neoplasm Malignant",
  "Astroblastoma",
  "Metastatic Glioma"
)

### 神经母细胞瘤相关PT####
histology_Neuroblastoma <- c(
  "Neuroblastoma",
  "Neuroblastoma Recurrent",
  "Neonatal Neuroblastoma",
  "Ganglioneuroblastoma",
  "Central Nervous System Neuroblastoma",
  "Aesthesioneuroblastoma"
)

### 淋巴瘤相关PT####
histology_Lymphoma <- c(
  # B细胞淋巴瘤
  "B-Cell Lymphoma",
  "B-Cell Lymphoma Refractory",
  "B-Cell Lymphoma Recurrent",
  "B-Cell Unclassifiable Lymphoma High Grade",
  "B-Cell Unclassifiable Lymphoma Low Grade",
  "B-Cell Lymphoma Unclassifiable",
  "High-Grade B-Cell Lymphoma",
  "High Grade B-Cell Lymphoma Refractory",
  "High Grade B-Cell Lymphoma Recurrent",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma Recurrent",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma Refractory",
  "Primary Mediastinal Large B-Cell Lymphoma",
  "Primary Mediastinal Large B-Cell Lymphoma Refractory",
  "Primary Mediastinal Large B-Cell Lymphoma Recurrent",
  "Primary Gastrointestinal Follicular Lymphoma",
  "B-Cell Small Lymphocytic Lymphoma",
  "B-Cell Small Lymphocytic Lymphoma Recurrent",
  "B-Cell Small Lymphocytic Lymphoma Refractory",
  # 套细胞淋巴瘤
  "Mantle Cell Lymphoma",
  "Mantle Cell Lymphoma Recurrent",
  "Mantle Cell Lymphoma Refractory",
  # 滤泡性淋巴瘤
  "Follicular Lymphoma",
  "Follicular Lymphoma Refractory",
  "Follicular Lymphoma Recurrent",
  "Follicle Centre Lymphoma, Follicular",
  "Follicle Centre Lymphoma Diffuse Small Cell Lymphoma",
  "Follicle Centre Lymphoma Diffuse Small Cell Lymphoma Recurrent",
  # T细胞淋巴瘤
  "T-Cell Lymphoma",
  "T-Cell Lymphoma Refractory",
  "T-Cell Lymphoma Recurrent",
  "Cutaneous T-Cell Lymphoma",
  "Cutaneous T-Cell Lymphoma Refractory",
  "Cutaneous T-Cell Lymphoma Recurrent",
  "Peripheral T-Cell Lymphoma Unspecified",
  "Peripheral T-Cell Lymphoma Unspecified Refractory",
  "Peripheral T-Cell Lymphoma Unspecified Recurrent",
  "Angioimmunoblastic T-Cell Lymphoma",
  "Angioimmunoblastic T-Cell Lymphoma Refractory",
  "Angioimmunoblastic T-Cell Lymphoma Recurrent",
  "Adult T-Cell Lymphoma/Leukaemia",
  "Adult T-Cell Lymphoma/Leukaemia Recurrent",
  "Adult T-Cell Lymphoma/Leukaemia Refractory",
  "Enteropathy-Associated T-Cell Lymphoma",
  "Intestinal T-Cell Lymphoma",
  "Intestinal T-Cell Lymphoma Recurrent",
  "Intestinal T-Cell Lymphoma Refractory",
  "Hepatosplenic T-Cell Lymphoma",
  "Angiocentric Lymphoma",
  "Angiocentric Lymphoma Recurrent",
  "Natural Killer-Cell Lymphoblastic Lymphoma",
  "Extranodal Nk/T-Cell Lymphoma, Nasal Type",
  # 霍奇金淋巴瘤
  "Hodgkin's Disease",
  "Hodgkin's Disease Nodular Sclerosis",
  "Hodgkin's Disease Refractory",
  "Hodgkin's Disease Recurrent",
  "Hodgkin's Disease Mixed Cellularity",
  "Hodgkin's Disease Mixed Cellularity Recurrent",
  "Hodgkin's Disease Mixed Cellularity Refractory",
  "Hodgkin's Disease Lymphocyte Predominance Type",
  "Hodgkin's Disease Lymphocyte Predominance",
  "Hodgkin's Disease Lymphocyte Depletion Type",
  "Hodgkin's Disease Lymphocyte Depletion",
  "Hodgkin's Disease Lymphocyte Depletion Type Refractory",
  "Hodgkin's Disease Nodular Sclerosis Recurrent",
  "Hodgkin's Disease Unclassifiable",
  "Nodular Lymphocyte Predominant Hodgkin Lymphoma",
  # 非霍奇金淋巴瘤
  "Non-Hodgkin's Lymphoma",
  "Non-Hodgkin's Lymphoma Recurrent",
  "Non-Hodgkin's Lymphoma Refractory",
  "Non-Hodgkin's Lymphoma Unspecified Histology Indolent",
  "Non-Hodgkin's Lymphoma Unspecified Histology Aggressive",
  "Non-Hodgkin's Lymphoma Unspecified Histology Aggressive Refractory",
  "Non-Hodgkin's Lymphoma Metastatic",
  "Non-Hodgkin's Lymphoma Transformed Recurrent",
  # 间变性大细胞淋巴瘤
  "Anaplastic Large-Cell Lymphoma",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types Refractory",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types Recurrent",
  "Breast Implant-Associated Anaplastic Large Cell Lymphoma",
  # 边缘区淋巴瘤
  "Marginal Zone Lymphoma",
  "Marginal Zone Lymphoma Recurrent",
  "Marginal Zone Lymphoma Refractory",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type)",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type) Recurrent",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type) Refractory",
  "Extranodal Marginal Zone B-Cell Lymphoma (Balt Type)",
  "Nodal Marginal Zone B-Cell Lymphoma",
  "Nodal Marginal Zone B-Cell Lymphoma Recurrent",
  "Nodal Marginal Zone B-Cell Lymphoma Refractory",
  "Splenic Marginal Zone Lymphoma",
  "Splenic Marginal Zone Lymphoma Recurrent",
  "Splenic Marginal Zone Lymphoma Refractory",
  # 伯基特淋巴瘤
  "Burkitt's Lymphoma",
  "Burkitt's Lymphoma Recurrent",
  "Burkitt's Lymphoma Refractory",
  # 其他淋巴瘤
  "Lymphoma",
  "Lymphoma Transformation",
  "Lymphoma Aids Related",
  "Gastrointestinal Lymphoma",
  "Ocular Lymphoma",
  "Lymphocytic Lymphoma",
  "Disseminated Large Cell Lymphoma",
  "Primary Effusion Lymphoma",
  "Primary Cardiac Lymphoma",
  "Double Hit Lymphoma",
  "Grey Zone Lymphoma",
  "Triple Hit Lymphoma",
  "Composite Lymphoma",
  "Immunoblastic Lymphoma",
  "Lymphoma Cutis",
  "Metastatic Lymphoma",
  "Leukaemic Lymphoma",
  "Epstein-Barr Virus Associated Lymphoma",
  "Malignant Lymphoma Unclassifiable High Grade",
  "Plasmablastic Lymphoma",
  "Central Nervous System Lymphoma",
  "Cutaneous Lymphoma",
  "Cutaneous B-Cell Lymphoma",
  "Primary Breast Lymphoma",
  "Thyroid B-Cell Lymphoma",
  # 华氏巨球蛋白血症
  "Lymphoplasmacytoid Lymphoma/Immunocytoma",
  "Lymphoplasmacytoid Lymphoma/Immunocytoma Recurrent",
  "Lymphoplasmacytoid Lymphoma/Immunocytoma Refractory",
  # 蕈样肉芽肿
  "Mycosis Fungoides",
  "Mycosis Fungoides Refractory"
)

### 白血病相关PT####
histology_Leukemia <- c(
  # 慢性淋巴细胞白血病
  "Chronic Lymphocytic Leukaemia",
  "Chronic Lymphocytic Leukaemia Recurrent",
  "Chronic Lymphocytic Leukaemia Refractory",
  "Chronic Lymphocytic Leukaemia (In Remission)",
  "Chronic Lymphocytic Leukaemia Transformation",
  "T-Cell Chronic Lymphocytic Leukaemia",
  # 慢性髓系白血病
  "Chronic Myeloid Leukaemia",
  "Chronic Myeloid Leukaemia (In Remission)",
  "Chronic Myeloid Leukaemia Recurrent",
  "Chronic Myeloid Leukaemia Transformation",
  "Philadelphia Positive Chronic Myeloid Leukaemia",
  "Blast Crisis In Myelogenous Leukaemia",
  # 急性淋巴细胞白血病
  "Acute Lymphocytic Leukaemia",
  "Acute Lymphocytic Leukaemia Recurrent",
  "Acute Lymphocytic Leukaemia Refractory",
  "Acute Lymphocytic Leukaemia (In Remission)",
  "Philadelphia Positive Acute Lymphocytic Leukaemia",
  "B Precursor Type Acute Leukaemia",
  "B-Cell Type Acute Leukaemia",
  "T-Cell Type Acute Leukaemia",
  "Precursor T-Lymphoblastic Leukaemia Acute",
  "Myeloblastoma",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia Refractory",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia Recurrent",
  "Precursor B-Lymphoblastic Lymphoma",
  "Precursor B-Lymphoblastic Lymphoma Recurrent",
  "Precursor B-Lymphoblastic Lymphoma Refractory",
  "Acute Biphenotypic Leukaemia",
  "Acute Bilineal Leukaemia",
  "Acute Undifferentiated Leukaemia",
  "Mature B-Cell Type Acute Leukaemia",
  "Burkitt's Leukaemia",
  # 其他白血病
  "Leukaemia",
  "Leukaemia Recurrent",
  "Leukaemia In Remission",
  "Lymphocytic Leukaemia",
  "Myeloid Leukaemia",
  "Myeloid Leukaemia In Remission",
  "Acute Leukaemia",
  "Acute Leukaemia In Remission",
  "Chronic Leukaemia",
  "Chronic Leukaemia In Remission",
  "Chronic Myelomonocytic Leukaemia",
  "Chronic Myelomonocytic Leukaemia (In Remission)",
  "Juvenile Chronic Myelomonocytic Leukaemia",
  "Hairy Cell Leukaemia",
  "Hairy Cell Leukaemia Recurrent",
  "T-Cell Prolymphocytic Leukaemia",
  "B-Cell Prolymphocytic Leukaemia",
  "Prolymphocytic Leukaemia",
  "Natural Killer-Cell Leukaemia",
  "Chronic Eosinophilic Leukaemia",
  "Eosinophilic Leukaemia",
  "Mastocytic Leukaemia",
  "Leukaemia Cutis",
  "Leukaemia Basophilic",
  "Leukaemia Granulocytic",
  "Leukaemia Monocytic",
  "Monocytic Leukaemia In Remission",
  "Lymphoid Leukaemia (In Remission)",
  "Leukaemia Recurrent Prophylaxis",
  "Lineage Switch Leukaemia",
  "Aleukaemic Leukaemia",
  "Central Nervous System Leukaemia",
  "Plasma Cell Leukaemia",
  "Plasma Cell Leukaemia In Remission"
)

### 浆细胞肿瘤相关PT####
histology_Plasma_Cell_Neoplasm <- c(
  "Plasma Cell Myeloma",
  "Plasma Cell Myeloma In Remission",
  "Plasma Cell Myeloma Recurrent",
  "Plasma Cell Myeloma Refractory",
  "Plasmacytoma",
  "Myeloma Cast Nephropathy"
)

### 骨髓增生异常/骨髓增殖性肿瘤相关PT####
histology_MDS_MPN <- c(
  "Myelodysplastic Syndrome",
  "Myeloproliferative Neoplasm",
  "Myelodysplastic Syndrome Unclassifiable",
  "Myelodysplastic Syndrome With Excess Blasts",
  "Myelodysplastic Syndrome With Multilineage Dysplasia",
  "Myelodysplastic Syndrome With Single Lineage Dysplasia",
  "Myelodysplastic Syndrome With Ringed Sideroblasts",
  "Myelodysplastic Syndrome Transformation",
  "Myeloproliferative Disorder"
)

### 生殖细胞肿瘤相关PT####
histology_Germ_Cell_Tumor <- c(
  "Germ Cell Cancer",
  "Germ Cell Cancer Metastatic",
  "Ovarian Germ Cell Cancer",
  "Ovarian Germ Cell Embryonal Carcinoma",
  "Ovarian Embryonal Carcinoma",
  "Extragonadal Primary Embryonal Carcinoma",
  "Testicular Cancer Metastatic",
  "Choriocarcinoma",
  "Choriocarcinoma Recurrent",
  "Metastatic Choriocarcinoma"
)

### 神经内分泌肿瘤相关PT####
histology_Neuroendocrine_Tumor <- c(
  "Neuroendocrine Carcinoma",
  "Neuroendocrine Carcinoma Metastatic",
  "Neuroendocrine Carcinoma Of The Skin",
  "Neuroendocrine Carcinoma Of The Oesophagus",
  "Mixed Adenoneuroendocrine Carcinoma",
  "Gastrointestinal Neuroendocrine Carcinoma",
  "Neuroendocrine Tumour Of The Lung Metastatic",
  "Pancreatic Neuroendocrine Tumour Metastatic",
  "Metastatic Carcinoid Tumour",
  "Endocrine Neoplasm Malignant",
  "Glucagonoma",
  "Metastatic Glucagonoma",
  "Malignant Neoplasm Of Islets Of Langerhans"
)

### 胚胎性肿瘤相关PT####
histology_Embryonal_Tumor <- c(
  "Medulloblastoma",
  "Medulloblastoma Recurrent",
  "Hepatoblastoma",
  "Hepatoblastoma Recurrent",
  "Pancreatoblastoma",
  "Nephroblastoma",
  "Retinoblastoma",
  "Congenital Retinoblastoma",
  "Extraocular Retinoblastoma",
  "Pinealoblastoma",
  "Primitive Neuroectodermal Tumour Metastatic"
)

### 间皮瘤相关PT####
histology_Mesothelioma <- c(
  "Malignant Neoplasm Of Pleura",
  "Malignant Neoplasm Of Pleura Metastatic"
)

### 胸腺瘤相关PT####
histology_Thymoma <- c(
  "Thymic Cancer Metastatic"
)

### 肝细胞癌相关PT####
histology_Hepatocellular_Carcinoma <- c(
  "Hepatic Cancer Metastatic"
  # 注意：Hepatoblastoma归类为Embryonal Tumor
)

### 肾细胞癌相关PT####
histology_Renal_Cell_Carcinoma <- c(
  "Renal Cell Carcinoma",
  "Metastatic Renal Cell Carcinoma",
  "Renal Cell Carcinoma Recurrent",
  "Non-Renal Cell Carcinoma Of Kidney",
  "Hereditary Leiomyomatosis Renal Cell Carcinoma",
  "Renal Cancer",
  "Renal Cancer Recurrent",
  "Renal Cancer Metastatic"
)

### 甲状腺癌相关PT####
histology_Thyroid_Carcinoma <- c(
  "Huerthle Cell Carcinoma",
  "Thyroid Cancer Metastatic"
  # 注意：Thyroid B-Cell Lymphoma归类为Lymphoma
)

### 基底细胞癌相关PT####
histology_Basal_Cell_Carcinoma <- c(
  "Basal Cell Carcinoma",
  "Basal Cell Carcinoma Metastatic",
  "Vulvar Basal Cell Carcinoma"
)

### 其他癌相关PT（能确定是癌但无法细分类型）####
histology_Other_Carcinoma <- c(
  # 大细胞肺癌
  "Large Cell Lung Cancer",
  "Large Cell Lung Cancer Recurrent",
  "Large Cell Lung Cancer Metastatic",
  # 类癌/唾液腺特殊类型
  "Mucoepidermoid Carcinoma",
  "Adenoid Cystic Carcinoma",
  "Adenoid Cystic Carcinoma Of Salivary Gland",
  "Acinic Cell Carcinoma Of Salivary Gland",
  "Mucoepidermoid Carcinoma Of Salivary Gland",
  "Carcinoma Ex-Pleomorphic Adenoma",
  "Myoepithelial Carcinoma Of Salivary Gland",
  # 肾上腺皮质癌
  "Adrenal Gland Cancer",
  "Adrenal Gland Cancer Metastatic",
  # 皮肤附属器癌
  "Sebaceous Carcinoma",
  "Eccrine Carcinoma",
  "Pilomatrix Carcinoma",
  "Trichoblastic Carcinoma",
  # 脉络丛癌
  "Choroid Plexus Carcinoma",
  # 成釉细胞癌
  "Ameloblastic Carcinoma",
  # 印戒细胞癌（无法确定原发灶）
  "Signet-Ring Cell Carcinoma",
  # NUT中线癌
  "Nut Midline Carcinoma"
)

### Unknown组织学类型（无法从PT名称确定组织学类型）####
histology_Unknown <- c(
  # 肺部 - 未细分类型
  "Non-Small Cell Lung Cancer",
  "Lung Neoplasm Malignant",
  "Lung Carcinoma Cell Type Unspecified",
  "Non-Small Cell Lung Cancer Recurrent",
  "Bronchial Carcinoma",
  "Non-Small Cell Lung Cancer Metastatic",
  "Lung Carcinoma Cell Type Unspecified Recurrent",
  "Lung Cancer Metastatic",
  "Her2 Mutant Non-Small Cell Lung Cancer",
  "Metastatic Bronchial Carcinoma",
  "Respiratory Tract Carcinoma In Situ",
  # 膀胱/泌尿道 - 非尿路上皮癌类型
  "Urethral Cancer",
  "Ureteric Cancer",
  "Ureteric Cancer Recurrent",
  "Ureteric Cancer Local",
  "Ureteric Cancer Metastatic",
  "Metastatic Carcinoma Of The Bladder",
  "Urethral Cancer Metastatic",
  "Malignant Urinary Tract Neoplasm Metastatic",
  # 肾脏 - 非肾细胞癌类型
  "Malignant Neoplasm Of Renal Pelvis",
  # 前列腺 - 泛称
  "Prostate Cancer Metastatic",
  "Malignant Neoplasm Of Seminal Vesicle",
  "Malignant Neoplasm Of Paraurethral Glands",
  # 宫颈
  "Cervix Cancer Metastatic",
  # 子宫
  "Endometrial Cancer",
  "Malignant Neoplasm Of Uterine Adnexa",
  "Uterine Carcinoma In Situ",
  "Endometrial Cancer Metastatic",
  "Metastatic Uterine Cancer",
  # 乳腺
  "Breast Cancer Metastatic",
  # 腹膜
  "Metastases To Peritoneum",
  "Peritoneal Cancer Index",
  # 卵巢 - 非生殖细胞类型
  "Fallopian Tube Cancer",
  "Ovarian Stromal Cancer",
  "Fallopian Tube Cancer Metastatic",
  "Ovarian Cancer Metastatic",
  # 肠道 - 泛称类型
  "Colorectal Cancer",
  "Rectosigmoid Cancer",
  "Anal Cancer",
  "Appendix Cancer",
  "Colorectal Cancer Metastatic",
  "Colorectal Cancer Recurrent",
  "Small Intestine Carcinoma",
  "Rectosigmoid Cancer Recurrent",
  "Rectosigmoid Cancer Metastatic",
  "Small Intestine Carcinoma Metastatic",
  "Small Intestine Carcinoma Recurrent",
  "Anal Cancer Recurrent",
  "Gastrointestinal Carcinoma",
  "Gastrointestinal Carcinoma In Situ",
  "Colorectal Carcinoma",
  "Gastrointestinal Cancer Metastatic",
  "Intestinal Metastasis",
  "Hereditary Non-Polyposis Colorectal Cancer Syndrome",
  "Her2 Positive Colorectal Cancer",
  "Recurrent N-Ras Mutation-Positive Colorectal Carcinoma",
  "Colon Cancer Metastatic",
  "Rectal Cancer Metastatic",
  "Appendix Cancer Metastatic",
  "Anal Cancer Metastatic",
  # 壶腹部
  "Malignant Neoplasm Of Ampulla Of Vater",
  "Malignant Neoplasm Papilla Of Vater",
  # 食管胃
  "Gastrooesophageal Cancer",
  "Gastrooesophageal Cancer Recurrent",
  "Oesophageal Cancer Metastatic",
  "Metastatic Gastric Cancer",
  # 胆道
  "Gallbladder Cancer",
  "Gallbladder Cancer Metastatic",
  "Gallbladder Cancer Recurrent",
  "Biliary Cancer Metastatic",
  # 胰腺 - 非神经内分泌类型
  "Pancreatic Carcinoma Metastatic",
  # 头颈部
  "Salivary Gland Cancer",
  "Throat Cancer",
  "Gingival Cancer",
  "Her2 Positive Salivary Gland Cancer",
  "Salivary Gland Cancer Recurrent",
  "Metastatic Salivary Gland Cancer",
  "Epiglottic Cancer",
  "Epiglottic Carcinoma",
  "Tracheal Cancer",
  "Postcricoid Cancer",
  "Ear Neoplasm Malignant",
  "External Ear Neoplasm Malignant",
  "Malignant Neoplasm Of Auricular Cartilage",
  "Head And Neck Cancer Metastatic",
  "Tongue Cancer Metastatic",
  "Sinus Cancer Metastatic",
  # 皮肤 - 泛称
  "Skin Cancer",
  "Carcinoma In Situ Of Skin",
  "Skin Cancer Metastatic",
  # 脑部 - 非胶质瘤类型
  "Haemangioblastoma",
  "Brain Cancer Metastatic",
  "Pineal Parenchymal Neoplasm Malignant",
  "Malignant Neoplasm Of Spinal Cord",
  "Pituitary Neoplasm Malignant Recurrent",
  "Pituitary Neoplasm Malignant",
  "Pituitary Cancer Metastatic",
  "Malignant Meningioma Metastatic",
  # 骨
  "Bone Cancer",
  "Bone Cancer Metastatic",
  "Chondroblastoma",
  # 软组织 - 泛称
  "Soft Tissue Cancer",
  "Retroperitoneal Cancer",
  "Retroperitoneal Neoplasm Metastatic",
  # CUP和泛称
  "Malignant Neoplasm Of Unknown Primary Site",
  "Neoplasm Malignant",
  "Carcinoma In Situ",
  "Malignant Neoplasm Progression",
  "Recurrent Cancer",
  "Refractory Cancer",
  "Cancer In Remission",
  "Targeted Cancer Therapy",
  "Cancer Pain",
  "Cancer Surgery",
  "Cancer Hormonal Therapy",
  "Cancer Screening",
  "Cancer Staging",
  "Cancer Fatigue",
  "Cancer With A High Tumour Mutational Burden",
  "Microsatellite Instability Cancer",
  "Ntrk Gene Fusion Cancer",
  "Pd-L1 Positive Cancer",
  "Mismatch Repair Cancer Syndrome",
  # Unknown Primary转移相关
  "Metastases To Bone",
  "Metastases To Lung",
  "Metastases To Liver",
  "Metastases To Central Nervous System",
  "Metastases To Lymph Nodes",
  "Metastases To Meninges",
  "Metastases To Spine",
  "Metastases To Kidney",
  "Metastases To The Mediastinum",
  "Metastases To Abdominal Cavity",
  "Metastases To Adrenals",
  "Metastases To Bone Marrow",
  "Metastases To Pleura",
  "Metastases To Nervous System",
  "Metastases To Eye",
  "Metastases To Skin",
  "Metastases To The Respiratory System",
  "Metastases To Breast",
  "Metastases To Gastrointestinal Tract",
  "Metastases To Biliary Tract",
  "Metastases To Bladder",
  "Metastases To Muscle",
  "Metastases To Neck",
  "Metastases To Ovary",
  "Metastases To Pancreas",
  "Metastases To Pelvis",
  "Metastases To Peripheral Nervous System",
  "Metastases To Pituitary Gland",
  "Metastases To Prostate",
  "Metastases To Retroperitoneum",
  "Metastases To Soft Tissue",
  "Metastases To Spleen",
  "Metastases To Stomach",
  "Metastases To Tonsils",
  "Metastases To Trachea",
  "Metastases To Uterus",
  "Metastases To Rectum",
  "Metastases To Spinal Cord",
  "Metastases To Thorax",
  "Metastases To Urinary Tract",
  "Metastases To Diaphragm",
  "Metastases To Heart",
  "Metastases To Large Intestine",
  "Metastases To Mouth",
  "Metastases To Oesophagus",
  "Metastases To Penis",
  "Metastases To Chest Wall",
  "Metastases To Thyroid",
  "Metastases To Reproductive Organ",
  "Metastases To Salivary Gland",
  "Metastases To Abdominal Wall",
  "Metastases To Gallbladder",
  "Metastases To Nasal Sinuses",
  "Metastases To Pharynx",
  "Metastases To Placenta",
  "Metastases To Small Intestine",
  "Metastases To Testicle",
  "Metastases To Vagina",
  "Metastasis",
  "Metastatic Neoplasm",
  "Metastatic Pain",
  "Metastatic Bone Disease Prophylaxis",
  "Infected Metastasis",
  "Hypercalcaemia Of Malignancy",
  "Anaemia Of Malignant Disease",
  "Polyneuropathy In Malignant Disease",
  "Second Primary Malignancy",
  "Maternal Cancer In Pregnancy",
  "Cancer Cells Urine Present",
  "Malignant Neoplasm Oligoprogression",
  "Peritoneal Carcinoma Metastatic",
  # 稀有肿瘤
  "Malignant Neoplasm Of Placenta",
  "Splenic Neoplasm Malignancy Unspecified",
  "Cardiac Neoplasm Malignant",
  "Malignant Neoplasm Of Thorax",
  "Abdominal Wall Neoplasm Malignant",
  "Ameloblastoma",
  "Anti-Melanoma Differentiation-Associated Protein 5 Antibody Positive",
  "Kaposi's Varicelliform Eruption",
  # 眼部 - 非黑色素瘤非视网膜母细胞瘤
  "Malignant Neoplasm Of Choroid",
  "Malignant Neoplasm Of Eye",
  "Malignant Neoplasm Of Conjunctiva",
  "Malignant Neoplasm Of Lacrimal Duct",
  "Malignant Neoplasm Of Lacrimal Gland",
  "Malignant Neoplasm Of Retina",
  "Malignant Neoplasm Of Orbit",
  "Malignant Neoplasm Of Eyelid",
  "Malignant Neoplasm Of Cornea",
  "Carcinoma In Situ Of Eye",
  "Ocular Cancer Metastatic",
  # 外阴阴道 - 泛称
  "Vulval Cancer",
  "Vulval Cancer Recurrent",
  "Vulval Cancer Metastatic",
  "Vaginal Cancer",
  "Vaginal Cancer Recurrent",
  "Genital Neoplasm Malignant Female",
  "Female Reproductive Tract Carcinoma In Situ",
  "Vaginal Cancer Metastatic",
  # 男性生殖器 - 泛称
  "Penile Cancer",
  "Penis Carcinoma",
  "Penis Carcinoma Metastatic",
  "Carcinoma In Situ Of Penis",
  "Genital Cancer Male",
  "Scrotal Cancer",
  "Epididymal Cancer",
  "Malignant Neoplasm Of Spermatic Cord",
  # 血液系统泛称
  "Haematological Malignancy"
)

# 良性/交界性肿瘤PT列表（用于排除）####
# 说明：这些PT应在cohort_mapping时被标记为exclude_benign
# 在2.1/2.2/2.3分析集处理脚本中的benign_pattern已更新以匹配这些PT
benign_tumor_pts <- c(
  # 明确良性肿瘤
  "Benign Mesothelioma",
  "Benign Pleural Mesothelioma",
  "Lipoblastoma",
  "Angiomyofibroblastoma",
  "Juvenile Melanoma Benign",
  # WHO Grade I / 交界性肿瘤（2026-01-18新增）
  "Haemangioblastoma",        # 血管母细胞瘤，WHO Grade I
  "Hemangioblastoma",         # 同上（美式拼写）
  "Ameloblastoma",            # 成釉细胞瘤，局部侵袭性良性
  "Chondroblastoma"           # 软骨母细胞瘤，良性骨肿瘤
)

# 非肿瘤疾病PT列表（被恶性模式误匹配）####
# 说明：这些PT包含"metastatic"但实际不是恶性肿瘤
non_malignant_pts <- c(
  "Metastatic Cutaneous Crohn's Disease",  # 克罗恩病皮肤表现
  "Calcification Metastatic",              # 转移性钙化（非肿瘤）
  # 2026-01-28 JADER最终映射：皮肌炎（自身免疫性疾病，非肿瘤）
  "Anti-Mda5-Dm Anti-Melanoma Differentiation-Associated Protein 5 Antibody-Positive Dermatomyositis",
  "Anti-Melanoma Differentiation-Associated Gene 5-Positive Dermatomyositis",
  "Clinically Amyopathic Dermatomyositis (Cadm",
  # 2026-01-28 JADER最终映射：检查方式（非疾病诊断）
  "Bronchoscopy"
)

# 第一层映射规则：Primary organs/tissues（基于first_pt）####
# 格式：organ_xxx <- c("Pt Name 1", "Pt Name 2", ...)
# 注意：first_pt已转为首字母大写格式（str_to_title）
# 说明：当 assigned_type 为 NA 或 "Others" 时，使用此规则

## Lung（肺）####
## 包含：NSCLC, SCLC, 肺癌相关
organ_Lung <- c(
  # 核心PT（2026-01-28添加）
  "Lung Adenocarcinoma",
  "Lung Squamous Cell Carcinoma",
  "Squamous Cell Carcinoma Of Lung",
  "Pulmonary Adenocarcinoma",
  "Non-Small Cell Lung Cancer",
  "Lung Neoplasm Malignant",
  "Small Cell Lung Cancer",
  "Lung Carcinoma Cell Type Unspecified",
  "Non-Small Cell Lung Cancer Recurrent",
  "Bronchial Carcinoma",
  "Small Cell Lung Cancer Extensive",
  "Non-Small Cell Lung Cancer Metastatic",
  "Lung Carcinoma Cell Type Unspecified Recurrent",
  "Large Cell Lung Cancer",
  "Sarcomatoid Carcinoma Of The Lung",
  "Lung Cancer Metastatic",
  "Her2 Mutant Non-Small Cell Lung Cancer",
  "Small Cell Lung Cancer Limited",
  "Bronchioloalveolar Carcinoma",
  "Large Cell Lung Cancer Recurrent",
  "Small Cell Lung Cancer Recurrent",
  "Small Cell Lung Cancer Metastatic",
  "Metastatic Bronchial Carcinoma",
  "Primary Pulmonary Melanoma",
  "Respiratory Tract Carcinoma In Situ",
  # 新增
  "Lung Adenocarcinoma Metastatic",
  "Lung Squamous Cell Carcinoma Metastatic",
  "Large Cell Lung Cancer Metastatic",
  "Neuroendocrine Tumour Of The Lung Metastatic",
  "Small Cell Lung Cancer Extensive Stage",
  # 2026-01-20 新增未映射PT
  "Lung Cancer",
  "Neuroendocrine Tumour Of The Lung",
  # 2026-01-28 JADER新增（ALK/ROS1/MET相关）
  "Alk-Positive Lung Cancer",
  "Alk Positive Lung Cancer",
  "Alk Lung Cancer",
  "Alk Fusion Gene-Positive Lung Cancer",
  "Ros1 Positive Lung Cancer",
  "Met Positive Lung Cancer",
  "Small Cell Lung Cancer",
  "Advanced Non-Small Cell Lung Cancer",
  "Alk-Positive Non-Small Cell Lung Cancer",
  "Alk Positive Non-Small Cell Lung Cancer",
  # 2026-01-28 JADER新增（长描述变体）
  "Unresectable, Advanced Or Recurrent Non-Small Cell Lung Cancer With Met Exon 14 Skipping Alterations",
  "Metex14 Skipping Alterations, Unresectable, Advanced Or Recurrent Non-Small Cell Lung Cancer",
  "Unresectable, Advanced Or Recurrent Non-Small Cell Lung Cancer",
  "Non-Small Cell Lung Cancer Associated With Interstitial Pneumonia",
  "Anaplastic Lymphoma Kinase Fusion Gene-Positive Unresectable Advanced And/Or Recurrent Non-Small Cell Lung Cancer With Resistance Or Intolerance To Alk Tyrosine Kinase Inhibitor(S",
  "Alk Fusion Gene-Positive Non-Small Cell Lung Cancer",
  "Alk Fusion Gene-Positive Unresectable Advanced And/Or Recurrent Non-Small Cell Lung Cancer",
  "Alk Fusion Gene-Positive Unresectable Advanced/Recurrent Non-Small Cell Lung Cancer Resistant Or Intolerant To Alk Tyrosine Kinase Inhibitors",
  "Lung Metastases",
  "Lung Metastasis",
  # 2026-01-28 JADER完整映射（肺癌相关）
  "Non-Small Cell Lung Cancer (Nsclc",
  "As A Preoperative Chemotherapy For Lung Cancer",
  "Non-Small Cell Cancer",
  "Non-Small Cell Lung Cancer With Her2 Mutation",
  "Non-Small Cell Lung Cancer With Met Exon 14 Skipping Alterations",
  "Akl Positive Lung Cancer",
  "Alk Fusion Gene Positive Lung Cancer",
  "Alk Fusion Gene-Positive Non-Small-Cell Lung Cancer (Nsclc",
  "Alk Fusion Gene-Positive Unresectable Advanced / Recurrent Non-Small Cell Lung Cancer Resistant Or Intolerant To Alk Tyrosine Kinase Inhibitors",
  "Alk Fusion Gene-Positive Unresectable Advanced And/Or Recurrent Non-Small Cell Lung Cancer With Resistance Or Intolerance To Alk Tyrosine Kinase Inhibitor(S",
  "Alk Fusion Gene-Positive Unresectable Advanced With Non-Small Cell Lung Cancer",
  "Alk Positib Lung Cancer",
  "Alk Positive Lung Cancer Treatment",
  "Alk Positive Lung Cancer/ Brain Metastases",
  "Egfr-Positive Lung Cancer",
  "Lung Carcinoma",
  "Left Lower Lung Lobe-Primary Small Cell Lung Cancer",
  "Left Superior Lobe Non-Small Cell Lung Cancer",
  "Lung Cancer In The Right Inferior Lobe",
  "Lung Cancer In The Right Lower Lobe Met Positive",
  "Lung Cancer Recurrence",
  "Lung Cancer Therapy",
  "Lung Cancer With Egfr Gene Mutation",
  "Lung Cancer With Met Exon 14 Skipping Mutation With Bone Metastases",
  "Lung Cancer With Met Positive",
  "Lung Cancer With Positive Met Exon 14 Skipping Alterations",
  "Malignant Neoplasm Of Right Upper Lobe (Alk Positive",
  "Met Exon 14-Skipping Mutation-Positive Non-Small Cell Lung Cancer",
  "Met Positive Nonsmall Cell Lung Cancer",
  "Metastatic Alk-Positive Pulmonary Adenocarcinoma",
  "Non Small Cell Lung Cancer",
  "Non Small Cell Lung Cancer With Met Exon 14 Skipping Mutation",
  "Non-Small Cell Lung Cancer (With Bone, Adrenal, And Other Metastases",
  "Non-Small Cell Lung Cancer Nos",
  "Non-Small Cell Lung Cancer With A Positive Alk Fusion Gene",
  "Non-Small Cell Lung Cancer With Met Gene Exon 14 Skipping Mutation",
  "Non-Small Cell Lung Cancer?",
  "Non-Small-Cell Lung Cancer",
  "Non-Small-Cell Lung Cancer ( Nsclc",
  "Postoperative Recurrent Ros1-Positive Lung Cancer",
  "Primary Pulmonary Squamous Cell Carcinoma",
  "Progressive Non-Small Cell Lung Cancer",
  "Pulmonary Sarcomatoid Carcinoma",
  "Right Lower Lobe Lung Cancer",
  "Right Lower Lobe Lung Cancer (Non-Small Cell Lung Cancer",
  "Right Lower Lobe Non-Small Cell Lung Cancer (",
  "Right Middle Lobe Lung Cancer (",
  "Right Pulmonary Hilum Lung Cancer",
  "Right Upper Lobe Lung Cancer",
  "Right Upper Lobe Lung Cancer Bone Metastasis",
  "Ros1 Fusion Gene Positive Non-Small Cell Lung Cancer",
  "Ros1 Fusion Gene-Positive Lung Cancer",
  "Ros1 Gene Positive Lung Cancer",
  "Ros1 Lung Cancer",
  "Ros1 Mutation Positive Lung Cancer",
  "Ros1-Positive Lung Cancer",
  "Ros1-Positive Lung Cancer Treatment",
  "Small Cell Lung Cancer (Sclc",
  "Treatment For Alk Lung Cancer",
  "Unresectable Advanced Or Recurrent Non-Small Cell Lung Cancer With Met Exon 14 Skipping Alterations.",
  "Unresectable Advanced Recurrent Non-Small Cell Lung Cancer With Met Exon 14 Skipping-Mutation",
  "Unresectable Advanced Recurrent Non-Small Cell Lung Cancer With Metastatic Exon 14 Skipping Alterations.",
  "Unresectable Advanced/Recurrent Non-Small Cell Lung Cancer",
  "Unresectable, Advanced Or Recurrent Non-Small Cell Lung Cancer With Met Exon 14 Skipping Alterations.",
  "Hypopharyngeal Cancer With Lung Metastasis",
  "Lung And Bone Metastasis",
  "Lung And Retroperitoneum Lymph Node Metastasis",
  "Lung Cancer (Small Intestinal Metastasis",
  "Lung Carcinoma (Pd-L1 95%, Ttf-1 Negative, With Lung And Adrenal Metastases",
  "Metastasis Of Lung Cancer",
  "Multiple Lung Metastases",
  "Denovo Lung Metastasis",
  "Large Intestinal Carcinoma Recurrent With Multiple Lung Tumors",
  "Treatment For Metastases To Bone (Lung Cancer",
  # 2026-01-28 JADER最终映射：肺癌原发的转移
  "Bone Lesion Due To Lung Cancer Metastases To Bone",
  "Bone Lesions Due To Lung Cancer Bone Metastases",
  "Her2 Mutation Positive Non-Small Cell Lung Cancer"
)

## Pleura（胸膜）####
## 包含：Mesothelioma, 胸膜癌
## 注意：良性间皮瘤已移至benign_tumor_pts
organ_Pleura <- c(
  # 核心PT（2026-01-28添加）
  "Pleural Mesothelioma Malignant",
  "Mesothelioma",
  "Mesothelioma Malignant",
  "Peritoneal Mesothelioma",
  "Epithelioid Mesothelioma",
  "Malignant Neoplasm Of Pleura",
  "Malignant Neoplasm Of Pleura Metastatic"
)

## Kidney（肾）####
## 包含：Renal cell carcinoma, 肾癌相关
organ_Kidney <- c(
  # 核心PT（2026-01-28添加）
  "Clear Cell Renal Cell Carcinoma",
  "Renal Cell Carcinoma",
  "Renal Cancer",
  "Kidney Cancer",
  "Metastatic Renal Cell Carcinoma",
  "Malignant Neoplasm Of Renal Pelvis",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter",
  "Renal Cell Carcinoma Recurrent",
  "Renal Cancer Recurrent",
  "Non-Renal Cell Carcinoma Of Kidney",
  "Transitional Cell Cancer Of Renal Pelvis And Ureter Metastatic",
  "Hereditary Leiomyomatosis Renal Cell Carcinoma",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Recurrent",
  "Renal Cancer Metastatic",
  "Nephroblastoma",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Localised",
  "Transitional Cell Cancer Of The Renal Pelvis And Ureter Regional",
  # 2026-01-20 新增未映射PT（不带The的版本）
  "Transitional Cell Cancer Of Renal Pelvis And Ureter",
  # 2026-01-28 JADER新增
  "Curatively Unresectable Or Metastatic Renal Cell Carcinoma",
  "Advanced Renal Cell Carcinoma",
  "Renal Cell Cancer",
  "Clear Cell Carcinoma",
  "Metastatic Renal Cancer",
  "Unresectable Renal Cancer",
  # 2026-01-28 JADER完整映射
  "Curatively Unresectable Or Metastatic Renal Cell Carcinoma.",
  "Radically Unresectable Or Metastatic Renal Cell Carcinoma",
  "Left Renal Cancer",
  "Right Renal Cancer",
  "Left Renal Cancer (Ct1bn0m1",
  "Right Renal Cancer T3b",
  "Right-Sided Renal Cancer",
  "Ct3an2m1 Renal Carcinoma",
  "Ct3anom1 Left Renal Cancer",
  "Ineradicable/Non-Resectable Metastatic Renal Cell Carcinoma",
  "Intractable Advanced Renal Cell Carcinoma",
  "Metastatic Renal Carcinoma",
  "Metastatic Renal Cell Carcinoma/ Curatively Unresectable Or Metastatic Renal Cell Carcinoma",
  "Metastatic Renal Cell Carcinoma/Left Renal Cancer Ct3an2m1",
  "Non-Clear Cell Renal Cell Carcinoma",
  "Poor Risk Advanced Renal Cell Carcinoma",
  "Progressive Renal Cancer",
  "Recurrence Of Renal Cell Carcinoma",
  "Recurrent Postoperative Renal Cancer (Poor Risk According To The Imdc",
  "Recurrent Renal Cancer",
  "Recurrent Renal Cell Carcinoma",
  "Renal  Carcinoma",
  "Renal Carcinoma",
  "Renal Cancer Progression",
  "Renal Cancer/Renal Cell Carcinoma",
  "Renal Cell Carcinoma (T3bn0m1,",
  "Renal Cell Carcinoma / Recurrence Of Left Para-Aortic Lymph Node Metastasis",
  "Renal Cell Carcinoma With Lung Metastases",
  "Renal Cell Carcinoma With Multiple Lung Metastases (Pt4n2m 1,",
  "Right Renal Cell Carcinoma",
  "Metastasis Renal Cancer",
  "Retroperitoneal Metastasis After Surgery Of Left Renal Cell Carcinoma",
  "Secondary Treatment Of Metastatic Renal Cell Carcinoma",
  "Translocation Renal Cell Carcinoma With Recurrence In The Left Retroperitoneal Soft Tissue",
  "Treatment Of Renal Cancer",
  "Unresectable Or Metastatic Renal Cell Carcinoma",
  "Xp11.2-Translocated Left Renal Cell Carcinoma",
  "Renal Pelvis Cancer",
  "Left Renal Pelvis Cancer",
  "Right Renal Pelvis Carcinoma",
  "Advanced Renal Pelvis Cancer",
  # 2026-01-28 JADER最终映射：双侧肾癌
  "Bilateral Renal Cancer"
)

## Bladder（膀胱）####
## 包含：Bladder/urinary tract cancer, 膀胱及泌尿道癌
organ_Bladder <- c(
  "Transitional Cell Carcinoma",
  "Transitional Cell Carcinoma Urethra",
  "Transitional Cell Carcinoma Recurrent",
  "Transitional Cell Carcinoma Metastatic",
  "Urethral Cancer",
  "Ureteric Cancer",
  "Urinary Tract Carcinoma In Situ",
  "Ureteric Cancer Recurrent",
  "Ureteric Cancer Local",
  "Ureteric Cancer Metastatic",
  # 新增
  "Bladder Transitional Cell Carcinoma Metastatic",
  "Metastatic Carcinoma Of The Bladder",
  "Urethral Cancer Metastatic",
  "Malignant Urinary Tract Neoplasm Metastatic",
  # 2026-01-20 新增未映射PT
  "Malignant Urinary Tract Neoplasm",
  # 2026-01-28 JADER完整映射
  "Upper Tract Urothelial Carcinoma",
  "Urothelial Carcinoma (Ureteric Cancer",
  "Urothelial Carcinoma In The Left Renal Pelvis",
  "High-Grade Urothelial Carcinoma (Uc",
  "Metastatic Right Ureteric Cancer",
  "Metastatic Ureteral Cancer (T4n2m1,",
  "Postoperative Recurrence Of Ureteric Cancer",
  "Right Ureteric Cancer",
  "Right Ureteric Cancer With Metastasis",
  "Unresectable Advanced Ureteral Cancer"
)

## Prostate（前列腺）####
## 包含：Prostate cancer
organ_Prostate <- c(
  "Prostate Cancer",  # 标准化后的名称
  "Prostate Cancer Metastatic",
  "Prostate Cancer Recurrent",
  "Prostate Cancer Refractory",
  "Malignant Neoplasm Of Seminal Vesicle",
  "Malignant Neoplasm Of Paraurethral Glands",
  # 2026-01-28 JADER最终映射：日文前列腺癌
  "Gs　4+4=8の前立腺癌、T4n0m1b、Adenocarcinoma、Crpc（去勢抵抗性前立腺がん）"
)

## Cervix（宫颈）####
## 包含：Cervical cancer
organ_Cervix <- c(
  "Cervix Cancer",  # 标准化后的名称
  "Cervix Cancer Metastatic",
  "Cervix Cancer Recurrent"
)

## Uterus（子宫）####
## 包含：Endometrial cancer, Uterine sarcoma, Choriocarcinoma
organ_Uterus <- c(
  "Endometrial Cancer",
  "Malignant Neoplasm Of Uterine Adnexa",
  "Uterine Carcinoma In Situ",
  # 新增
  "Endometrial Cancer Metastatic",
  "Metastatic Uterine Cancer",
  "Metastatic Choriocarcinoma",
  "Choriocarcinoma",
  "Choriocarcinoma Recurrent"
)

## Breast（乳腺）####
## 包含：Breast cancer
organ_Breast <- c(
  "Breast Cancer",  # 标准化后的名称
  "Breast Cancer Metastatic",
  "Breast Cancer Recurrent",
  "Primary Breast Lymphoma",
  "Breast Implant-Associated Anaplastic Large Cell Lymphoma",
  # 2026-01-28 JADER完整映射
  "Invasive Ductal Carcinoma Iii B",
  "Invasive Lobular Carcinoma",
  "Right Invasive Ductal Carcinoma Her2 2+",
  "再発乳癌(Invasive Ductal Carcinoma Scirrhous Type"
)


## Ovary（卵巢）####
## 包含：Ovarian cancer, 输卵管癌
organ_Ovary <- c(
  "Ovarian Cancer",  # 标准化后的名称
  "Ovarian Cancer Metastatic",
  "Ovarian Cancer Recurrent",
  "Fallopian Tube Cancer",
  "Fallopian Tube Cancer Metastatic",
  "Ovarian Germ Cell Cancer",
  "Ovarian Germ Cell Embryonal Carcinoma",
  "Ovarian Stromal Cancer",
  "Ovarian Embryonal Carcinoma",
  # 2026-01-28 JADER完整映射
  "Primary Atypical Serous Carcinoma Of The Right Fallopian Tube",
  "Multiple Bone, Liver, And Ovarian Metastases",
  "Peritoneal Sereous Carcinoma"
)

## Bowel（肠道）####
## 包含：Colorectal cancer, Small bowel cancer, Appendiceal cancer, Anal cancer
organ_Bowel <- c(
  "Colorectal Cancer",
  "Rectosigmoid Cancer",
  "Colorectal Adenocarcinoma",
  "Anal Cancer",
  "Appendix Cancer",
  "Anal Squamous Cell Carcinoma",
  "Colorectal Cancer Metastatic",
  "Colorectal Cancer Recurrent",
  "Small Intestine Carcinoma",
  "Small Intestine Adenocarcinoma",
  "Adenocarcinoma Of Appendix",
  "Intestinal Adenocarcinoma",
  "Rectosigmoid Cancer Recurrent",
  "Rectosigmoid Cancer Metastatic",
  "Small Intestine Carcinoma Metastatic",
  "Small Intestine Carcinoma Recurrent",
  "Anal Cancer Recurrent",
  "Gastrointestinal Carcinoma",
  "Gastrointestinal Adenocarcinoma",
  "Gastrointestinal Carcinoma In Situ",
  "Colorectal Carcinoma",
  "Gastrointestinal Cancer Metastatic",
  "Mucinous Adenocarcinoma Of Appendix",
  "Intestinal Metastasis",
  "Hereditary Non-Polyposis Colorectal Cancer Syndrome",
  "Her2 Positive Colorectal Cancer",
  "Recurrent N-Ras Mutation-Positive Colorectal Carcinoma",
  # 新增（带状态修饰词的变体）
  "Colon Cancer Metastatic",
  "Rectal Cancer Metastatic",
  "Appendix Cancer Metastatic",
  "Anal Cancer Metastatic",
  # 2026-01-28 JADER新增（日英转换后的英文名称）
  "Large Intestine Carcinoma",
  "Caecal Cancer",
  "Metastatic Colorectal Cancer",
  "Cecal Cancer",
  "Large Intestine Cancer",
  "Brafv600 Mutant Colorectal Cancer",
  "Rectosigmoid Cancer",
  "Sigmoid Cancer",
  "Sigmoid Carcinoma",
  "Unresectable, Progressive And Recurrent Large Intestine Carcinoma",
  "Peritoneal Cancer",
  # 2026-01-28 JADER完整映射
  "Advanced And Recurrent Large Intestine Carcinoma",
  "Anal Canal Cancer",
  "Anal Squamous Cell Carcinoma (With Lymph Node Metastases",
  "Appendiceal Carcinoma With Signet-Ring Cell Carcinoma",
  "Colorectal Cancer (Braf Mutation",
  "Colorectal Cancer (Sigmoid Colon Cancer",
  "Colorectal Cancer (With Lymph Node Metastasis",
  "Large Intestine Carcinoma Braf Mutation Type",
  "Large Intestine Carcinoma Non-Resectable",
  "Large Intestine Carcinoma With Metastatic Pulmonary Tumour",
  "Large Intestine Carcinoma/ Colorectal Cancer",
  "Metachronous Multiple Large Intestine Carcinoma",
  "Metastatic Colorectal Cancer (Mcrc",
  "Metastatic Colorectal Carcinoma",
  "Metastatic Large Intestine Carcinoma",
  "Metastatic Large Intestine Carcinoma With Braf 600e Mutation Positive",
  "Mutant Colorectal Metastasis",
  "Non-Resectable Large Intestine Carcinoma With Braf V600e Mutation",
  "Obstructive Colorectal Cancer",
  "Primary Colorectal Cancer With Distant Metastases",
  "Primary Ileal Cancer (Pt3n1am0, P",
  "Progressed After Chemotherapy, Unresectable, Advanced/Recurrent Colorectal Cancer With Msi-High",
  "Progressive And Recurrent Large Intestine Carcinoma",
  "Progressive Unresectable Cecal Cancer",
  "Ras Wild-Type Colorectal Cancer",
  "Recurrent Colorectal Cancer",
  "Recurrent Large Intestine Carcinoma",
  "Resectable Braf V600e Mutant Colorectal Metastasis",
  "Sigmoid Cancer Metastatic",
  "Sigmoid Colorectal Cancer",
  "Squamous Cell Carcinoma Of The Anal Canal",
  "Squamous Cell Carcinoma Of The Rectosigmoid Junction With Synchronous Metastases",
  "Metastatic Squamous Cell Carcinoma Of The Anal Canal",
  "Unresectable Advanced Or Recurrent Colorectal Cancer",
  "Unresectable Large Intestine Carcinoma",
  "Unresectable Recurrent Or Metastatic Large Intestine Carcinoma (Splenic Flexure",
  "Villous Cancer",
  # 2026-01-28 JADER最终映射（BRAF变体）
  "Braf Mutant Colorectal Cancer",
  "Braf V600E Mutant Unresectable Advanced/Recurrent Colorectal Cancer",
  "Braf V600E Mutated Colorectal Cancer",
  "Braf V600E Mutation-Positive Unresectable/Metastatic Colorectal Cancer",
  "Braf V600E Positive Colorectal Cancer",
  "Colon Adenocarcinoma",
  # 2026-01-28 JADER最终映射（大肠癌BRAF变体补充）
  "Braf Gene Mutation Positive Large Intestine Carcinoma",
  "Braf Mutant Large Intestine Carcinoma",
  "Braf Mutation-Positive Metastatic Colorectal Cancer",
  "Braf V600e Mutant Metastatic Colorectal Cancer",
  "Braf V600e Mutation Positive Unresectable, Advanced Or Recurrent Large Intestine Carcinoma",
  "Braf V600e Mutation-Positive Non-Resectable Large Intestine Carcinoma",
  "Braf V600e-Mutant Metastatic Colorectal Cancer",
  "Brafv600 Mutant Metastatic Colorectal Cancer",
  "Brafv600e Mutation Non-Resectable Large Intestine Carcinoma"
)


## Esophagus_stomach（食管/胃）####
## 包含：食管胃癌
organ_Esophagus_stomach <- c(
  # 核心PT（2026-01-28添加）
  "Gastric Cancer",
  "Gastric Carcinoma",
  "Gastric Adenocarcinoma",
  "Adenocarcinoma Gastric",
  "Gastric Cancer Metastatic",
  "Gastric Cancer Recurrent",
  "Stomach Cancer",
  "Stomach Carcinoma",
  "Her2 Positive Gastric Cancer",
  "Oesophageal Cancer",
  "Oesophageal Carcinoma",
  "Esophageal Cancer",
  "Esophageal Carcinoma",
  "Oesophageal Squamous Cell Carcinoma",
  "Esophageal Squamous Cell Carcinoma",
  "Oesophageal Adenocarcinoma",
  "Esophageal Adenocarcinoma",
  "Gastrooesophageal Cancer",
  "Gastrooesophageal Cancer Recurrent",
  "Neuroendocrine Carcinoma Of The Oesophagus",
  # 新增
  "Oesophageal Cancer Metastatic",
  "Oesophageal Squamous Cell Carcinoma Metastatic",
  "Metastatic Gastric Cancer",
  # 2026-01-28 JADER新增
  "Advance Gastro-Oesophageal Cancer",
  "Advanced Colorectal Cancer (Crc",
  # 2026-01-28 JADER完整映射
  "Advance Gastro-Esophageal Cancer",
  "Advance Gastro-Oesophageal Cancer (",
  "Gastroesophageal Junction Cancer (With Metastasis",
  "Gastrooesophageal Junction Cancer",
  "Lymph Node Carcinoma Of The Thoracic Esophagus.",
  "T4 Lymph Node Carcinoma Of The Thoracic Esophagus",
  "Unresectable Advanced Gastroesophageal Junction Cancer (With Paraaortic Lymph Node Metastasis",
  "胃体部癌(Adenocarcinoma",
  "胃癌",
  "胸部食道腺癌(Lt,半周, Type3, Adenocarcinoma, Ct2n2m0,"
)

## Liver（肝）####
## 包含：Liver cancer, 肝母细胞瘤
organ_Liver <- c(
  # 核心PT（2026-01-28添加）
  "Hepatocellular Carcinoma",
  "Hepatocellular Carcinoma Recurrent",
  "Hepatocellular Carcinoma Metastatic",
  "Hepatic Cancer",
  "Liver Cancer",
  "Liver Cancer Recurrent",
  "Liver Cancer Metastatic",
  "Hepatoma",
  "Hepatoblastoma",
  "Hepatoblastoma Recurrent",
  # 新增
  "Hepatic Cancer Metastatic",
  # 2026-01-28 JADER完整映射
  "Liver Metastasis",
  "Multiple Liver Metastases"
)

## Biliary_tract（胆道）####
## 包含：Biliary tract cancer, 胆囊癌
organ_Biliary_tract <- c(
  "Gallbladder Cancer",
  "Gallbladder Adenocarcinoma",
  "Gallbladder Cancer Metastatic",
  "Gallbladder Cancer Recurrent",
  "Gallbladder Squamous Cell Carcinoma",
  # 新增
  "Biliary Cancer Metastatic",
  # 2026-01-28 JADER完整映射
  "Adenocarcinoma Of The Gallbladder",
  "Gallbladder Adenosquamous Carcinoma",
  "Locally Advanced Gallbladder Cancer",
  "Locally Advanced Neuroendocrine Carcinoma Of The Gallbladder",
  "Neuroendocrine Carcinoma Of Gallbladder"
)

## Pancreas（胰腺）####
## 包含：Pancreatic cancer, Glucagonoma
organ_Pancreas <- c(
  # 核心PT（2026-01-28添加）
  "Pancreatic Carcinoma",
  "Pancreatic Cancer",
  "Pancreatic Adenocarcinoma",
  "Pancreatic Carcinoma Recurrent",
  "Pancreatoblastoma",
  "Malignant Neoplasm Of Islets Of Langerhans",
  "Pancreatic Neuroendocrine Tumour Metastatic",
  # 新增
  "Pancreatic Carcinoma Metastatic",
  "Metastatic Glucagonoma",
  "Glucagonoma",
  # 2026-01-20 新增未映射PT
  "Pancreatic Neuroendocrine Tumour",
  # 2026-01-28 JADER完整映射
  "Metastatic Pancreatic Insulinoma"
)

## Head_neck（头颈部）####
## 包含：Salivary gland cancer, Head/neck cancer (other)
organ_Head_neck <- c(
  "Salivary Gland Cancer",
  "Throat Cancer",
  "Gingival Cancer",
  "Her2 Positive Salivary Gland Cancer",
  "Mucoepidermoid Carcinoma",
  "Adenoid Cystic Carcinoma",
  "Adenoid Cystic Carcinoma Of Salivary Gland",
  "Adenocarcinoma Of Salivary Gland",
  "Squamous Cell Carcinoma Of The Parotid Gland",
  "Acinic Cell Carcinoma Of Salivary Gland",
  "Salivary Gland Cancer Recurrent",
  "Metastatic Salivary Gland Cancer",
  "Mucoepidermoid Carcinoma Of Salivary Gland",
  "Epiglottic Cancer",
  "Epiglottic Carcinoma",
  "Carcinoma Ex-Pleomorphic Adenoma",
  "Myoepithelial Carcinoma Of Salivary Gland",
  "Tracheal Cancer",
  "Postcricoid Cancer",
  "Ear Neoplasm Malignant",
  "External Ear Neoplasm Malignant",
  "Malignant Neoplasm Of Auricular Cartilage",
  # 新增
  "Head And Neck Cancer Metastatic",
  "Tongue Cancer Metastatic",
  "Sinus Cancer Metastatic",
  # 2026-01-20 新增未映射PT
  "Sinus Cancer",
  # 2026-01-28 JADER新增
  "Mouth Cancer",
  "Mouth Floor Cancer",
  "Maxillary Gingival Cancer",
  "Cancer Of Maxillary Antrum",
  "Recurrent Mouth Cancer",
  "Recurrent Mouth Floor Cancer",
  "Buccal Mucosal Squamous Cell Carcinoma",
  "Mandibular Gingival Cancer",
  "Supraglottic Cancer",
  # 2026-01-28 JADER完整映射
  "Advanced Sphenoid Sinus Cancer",
  "Anaplastic Adenocarcinoma At The Submandibular Gland With Systemic Metastases",
  "Carcinoma Of The External Auditory Canal",
  "Carcinoma Of The Left Mandibular Gingiva",
  "Gingival Cancer Of The Right Mandible",
  "Left Lower Gum Cancer",
  "Left Mandibular Gingival Cancer",
  "Left Mandibular Gingival Squamous Cell Carcinoma",
  "Left Maxillary Cancer",
  "Left Maxillary Gingival Cancer",
  "Left Mouth Floor Cancer (T4an3bm0",
  "Left Parotid Cancer",
  "Locally Advanced Maxillary Sinus Cancer (Squamous Cell Carcinoma",
  "Locally Recurrent Lower Gingival Cancer",
  "Lower Gingival Squamous Cell Carcinoma",
  "Maxillary Carcinoma",
  "Metastatic Left Maxillary Gingival Cancer",
  "Mouth Cancer (Floor Of Mouth",
  "Mouth Cancer (Right Lower Gingiva",
  "Non Resectable And Locally Advanced Squamous Cell Carcinoma Of Right Lower Gum",
  "Parotid Gland Cancer",
  "Primary Squamous Cell Carcinoma Of The Left Buccal Mucosa",
  "Recurrent And Metastatic Mouth Cancer.",
  "Recurrent Buccal Mucosa Cancer",
  "Recurrent Cancer Of The Right Floor Of Mouth (Well Differentiated Squamous Cell Carcinoma",
  "Recurrent Left Maxillary Gingival Cancer",
  "Recurrent Salivary Gland Carcinoma",
  "Right Mandibular Gingiva Cancer",
  "Right Maxillary Gingival Cancer",
  "Right Mouth Floor Cancer",
  "Right Palatal Squamous Cell Carcinoma",
  "Right Parotid Gland Cancer",
  "Right Submandibular Gland Carcinoma",
  "Squamous Cell Carcinoma At The Left Lower Gum",
  "Squamous Cell Carcinoma At The Mouth Floor",
  "Squamous Cell Carcinoma In The Left Mouth Floor",
  "Squamous Cell Carcinoma Of Left Lower Gum",
  "Squamous Cell Carcinoma Of The External Auditory Canal",
  "Squamous Cell Carcinoma Of The Maxillary Sinus",
  "Submandibular Gland Cancer",
  "Tongue Cancer (Submandibular Recurrence",
  "Tongue Cancer With Lung Metastases",
  "Tracheal Squamous Cell Carcinoma",
  # 2026-01-28 JADER最终映射：口腔黏膜癌
  "Buccal Mucous Membrane Squamous Cell Carcinoma"
)

## Skin（皮肤）####
## 包含：Melanoma, Basal cell carcinoma, Skin cancer (other)
organ_Skin <- c(
  # 核心PT（2026-01-28添加）
  "Malignant Melanoma",
  "Melanoma",
  "Melanoma Recurrent",
  "Melanoma Metastatic",
  "Cutaneous Melanoma",
  "Basal Cell Carcinoma",
  "Skin Cancer",
  "Squamous Cell Carcinoma Of Skin",
  "Neuroendocrine Carcinoma Of The Skin",
  "Skin Squamous Cell Carcinoma Recurrent",
  "Skin Squamous Cell Carcinoma Metastatic",
  "Basal Cell Carcinoma Metastatic",
  "Cutaneous Lymphoma",
  "Cutaneous B-Cell Lymphoma",
  "Carcinoma In Situ Of Skin",
  "Basosquamous Carcinoma",
  "Basosquamous Carcinoma Of Skin",
  "Sebaceous Carcinoma",
  "Eccrine Carcinoma",
  "Pilomatrix Carcinoma",
  "Trichoblastic Carcinoma",
  # 新增
  "Skin Cancer Metastatic",
  "Metastatic Malignant Melanoma",
  "Urethral Melanoma Metastatic",
  # 2026-01-20 新增未映射PT
  "Skin Squamous Cell Carcinoma",
  "Urethral Melanoma",
  # 2026-01-28 JADER新增
  "Merkel Cell Carcinoma",
  "Curatively Unresectable Merkel Cell Carcinoma",
  "Metastatic Merkel Cell Carcinoma",
  # 2026-01-28 JADER完整映射
  "Curatively Unresectable Metastatic Merkel Cell Carcinoma",
  "Merkel Cell Carcinoma Of Right Lower Leg",
  "Merkel Cell Carcinoma Of The Upper Eyelid",
  "Metastatic Combined Merkel Cell Carcinoma",
  "Unresectable Merkel Cell Carcinoma"
)

## Brain（脑）####
## 包含：Glioma, 脑肿瘤相关
organ_Brain <- c(
  # 核心PT（2026-01-28添加）
  "Glioblastoma",
  "Glioblastoma Multiforme",
  "Glioma",
  "Glioma Malignant",
  "Astrocytoma",
  "Astrocytoma Malignant",
  "Brain Neoplasm Malignant",
  "Medulloblastoma",
  "Medulloblastoma Recurrent",
  "Central Nervous System Lymphoma",
  "Aesthesioneuroblastoma",
  "Haemangioblastoma",
  "Astroblastoma",
  "Brain Cancer Metastatic",
  "Central Nervous System Melanoma",
  "Central Nervous System Leukaemia",
  "Choroid Plexus Carcinoma",
  "Pinealoblastoma",
  "Pineal Parenchymal Neoplasm Malignant",
  "Malignant Neoplasm Of Spinal Cord",
  "Central Nervous System Neuroblastoma",
  "Pituitary Neoplasm Malignant Recurrent",
  # 新增
  "Pituitary Neoplasm Malignant",
  "Pituitary Cancer Metastatic",
  "Metastatic Glioma",
  "Primitive Neuroectodermal Tumour Metastatic",
  "Malignant Meningioma Metastatic",
  # 2026-01-20 新增未映射PT
  "Brain Cancer",
  "Pituitary Cancer",
  "Primitive Neuroectodermal Tumour",
  "Malignant Meningioma",
  # 2026-01-28 JADER完整映射
  "Anaplastic Astroblastoma",
  "Medulloblastoma/Posterior Fossa Primitive Neuroectodermal Tumor",
  "Relapsed Medulloblastoma"
)

## PNS（周围神经系统）####
## 包含：Neuroblastoma, 神经内分泌肿瘤
organ_PNS <- c(
  "Neuroblastoma",
  "Neuroblastoma Recurrent",
  "Neuroendocrine Carcinoma",
  "Ganglioneuroblastoma",
  "Neuroendocrine Carcinoma Metastatic",
  "Mixed Adenoneuroendocrine Carcinoma",
  "Gastrointestinal Neuroendocrine Carcinoma",
  "Neuroendocrine Tumour Of The Lung Metastatic",
  "Neonatal Neuroblastoma",
  "Endocrine Neoplasm Malignant",
  # 2026-01-28 JADER完整映射
  "Colorectal Neuroendocrine Carcinoma",
  "Large Cell Neuroendocrine Carcinoma",
  "Large-Cell Neuroendocrine Carcinoma",
  "Metastatic Neuroendocrine Carcinoma",
  "Neuroendocrine Carcinoma (Nec"
)

## Bone（骨）####
## 包含：Bone cancer, Chondrosarcoma
organ_Bone <- c(
  "Bone Cancer",
  "Osteosarcoma Metastatic",
  "Bone Cancer Metastatic",
  "Chondroblastoma",
  # 新增
  "Osteosarcoma Recurrent",
  "Chondrosarcoma Metastatic",
  # 2026-01-20 新增未映射PT
  "Osteosarcoma",
  "Chondrosarcoma",
  # 2026-01-28 JADER新增
  "Metastatic Osteosarcoma",
  "Maxillary Cancer",
  # 2026-01-28 JADER完整映射
  "Bone Metastases",
  "Metastatic Bone Tumour",
  "Metastatic Bone Tumors",
  "Multiple Bone Metastases",
  "Multiple Bone Metastasis",
  "Pulmonary Metastasis Of Osteosarcoma",
  "Sternum Metastas",
  "Systemic Bone Metastasis"
)

## Soft_tissue（软组织）####
## 包含：Soft tissue sarcoma, GIST, Kaposi's Sarcoma, Ewing's Sarcoma
organ_Soft_tissue <- c(
  "Ewing's Sarcoma",
  "Kaposi's Sarcoma",
  "Gastrointestinal Stromal Cancer",
  "Ewing's Sarcoma Recurrent",
  "Kaposi's Sarcoma Aids Related",
  "Ewing-Like Sarcoma",
  "Extra-Osseous Ewing's Sarcoma",
  "Kaposi's Sarcoma Classical Type",
  "Ewing's Sarcoma Metastatic",
  "Kaposi Sarcoma Inflammatory Cytokine Syndrome",
  "Immune Reconstitution Inflammatory Syndrome Associated Kaposi's Sarcoma",
  "Soft Tissue Cancer",
  "Retroperitoneal Cancer",
  # 新增
  "Sarcoma Metastatic",
  "Liposarcoma Metastatic",
  "Leiomyosarcoma Metastatic",
  "Extra-Osseous Ewing's Sarcoma Metastatic",
  "Angiosarcoma Metastatic",
  "Retroperitoneal Neoplasm Metastatic",
  "Malignant Fibrous Histiocytoma Metastatic",
  "Malignant Haemangiopericytoma Metastatic",
  "Extraskeletal Chondrosarcoma Metastatic",
  "Extraskeletal Osteosarcoma Metastatic",
  "Malignant Mesenchymoma Metastatic",
  # 2026-01-20 新增未映射PT
  "Retroperitoneal Neoplasm",
  "Malignant Fibrous Histiocytoma",
  "Sarcoma",
  "Malignant Haemangiopericytoma",
  "Malignant Mesenchymoma",
  "Extraskeletal Chondrosarcoma",
  # 2026-01-28 JADER新增
  "Gastrointestinal Stromal Sarcoma",
  "Extraskeletal Osteosarcoma",
  # 2026-01-28 JADER完整映射
  "Ewing Sarcoma/Bcor-Ccnb3 Sarcoma Arising In The Proximal Tibia",
  "Gastrointestinal Stromal Tumor",
  "Gastrointestinal Stromal Tumor Metastatic",
  "Intimal Sarcoma With Osteosarcomatous Differentiation",
  "Metastatic Esophageal Gist",
  "Metastatic Gastrointestinal Stromal Tumors",
  "Metastatic Solitary Fibrous Tumor"
)

## Thymus（胸腺）####
## 包含：Thymic tumor
organ_Thymus <- c(
  "Thymic Cancer Metastatic",
  # 2026-01-28 JADER完整映射
  "Thymic Mucinous Adenocarcinoma"
)

## Thyroid（甲状腺）####
## 包含：Thyroid cancer
organ_Thyroid <- c(
  "Huerthle Cell Carcinoma",
  "Thyroid B-Cell Lymphoma",
  # 新增
  "Thyroid Cancer Metastatic"
)

## CUP（原发灶不明癌）####
## 包含：Cancer of unknown primary, 不明原发癌，转移到某部位（无法判断原发灶）
## 注：合并了原CUP和Unknown_Primary两个类别
organ_CUP <- c(
  # 原CUP类别
  "Malignant Neoplasm Of Unknown Primary Site",
  "Neoplasm Malignant",
  "Adenocarcinoma",
  "Squamous Cell Carcinoma",
  # 2026-01-20 新增未映射PT（泛化描述，无法确定器官）
  "Gastrointestinal Cancer",
  "Cancer",
  "Adenosquamous Cell Carcinoma",
  "Sarcomatoid Carcinoma",
  "Small Cell Carcinoma",
  "Signet-Ring Cell Carcinoma",
  "Carcinoma In Situ",
  "Adenocarcinoma Metastatic",
  "Metastatic Squamous Cell Carcinoma",
  "Nut Midline Carcinoma",
  "Malignant Neoplasm Progression",
  "Recurrent Cancer",
  "Refractory Cancer",
  "Cancer In Remission",
  "Targeted Cancer Therapy",
  "Cancer Pain",
  "Cancer Surgery",
  "Cancer Hormonal Therapy",
  "Cancer Screening",
  "Cancer Staging",
  "Cancer Fatigue",
  "Cancer With A High Tumour Mutational Burden",
  "Microsatellite Instability Cancer",
  "Ntrk Gene Fusion Cancer",
  "Pd-L1 Positive Cancer",
  "Mismatch Repair Cancer Syndrome",
  # 原Unknown_Primary类别（转移到各器官/部位）
  "Metastases To Bone",
  "Metastases To Lung",
  "Metastases To Liver",
  "Metastases To Central Nervous System",
  "Metastases To Lymph Nodes",
  "Metastases To Meninges",
  "Metastases To Spine",
  "Metastases To Kidney",
  "Metastases To The Mediastinum",
  "Metastases To Abdominal Cavity",
  "Metastases To Adrenals",
  "Metastases To Bone Marrow",
  "Metastases To Pleura",
  "Metastases To Nervous System",
  "Metastases To Eye",
  "Metastases To Skin",
  "Metastases To The Respiratory System",
  "Metastases To Breast",
  "Metastases To Gastrointestinal Tract",
  "Metastases To Biliary Tract",
  "Metastases To Bladder",
  "Metastases To Muscle",
  "Metastases To Neck",
  "Metastases To Ovary",
  "Metastases To Pancreas",
  "Metastases To Pelvis",
  "Metastases To Peripheral Nervous System",
  "Metastases To Pituitary Gland",
  "Metastases To Prostate",
  "Metastases To Retroperitoneum",
  "Metastases To Soft Tissue",
  "Metastases To Spleen",
  "Metastases To Stomach",
  "Metastases To Tonsils",
  "Metastases To Trachea",
  "Metastases To Uterus",
  "Metastases To Peritoneum",
  "Metastases To Rectum",
  "Metastases To Spinal Cord",
  "Metastases To Thorax",
  "Metastases To Urinary Tract",
  "Metastases To Diaphragm",
  "Metastases To Heart",
  "Metastases To Large Intestine",
  "Metastases To Mouth",
  "Metastases To Oesophagus",
  "Metastases To Penis",
  "Metastases To Chest Wall",
  "Metastases To Thyroid",
  "Metastases To Reproductive Organ",
  "Metastases To Salivary Gland",
  "Metastases To Abdominal Wall",
  "Metastases To Gallbladder",
  "Metastases To Nasal Sinuses",
  "Metastases To Pharynx",
  "Metastases To Placenta",
  "Metastases To Small Intestine",
  "Metastases To Testicle",
  "Metastases To Vagina",
  # 转移相关泛称
  "Metastasis",
  "Metastatic Neoplasm",
  "Metastatic Carcinoid Tumour",
  "Metastatic Pain",
  "Metastatic Bone Disease Prophylaxis",
  "Metastatic Lymphoma",
  "Infected Metastasis",
  # 恶性肿瘤相关症状/并发症（无法确定原发灶）
  "Hypercalcaemia Of Malignancy",
  "Anaemia Of Malignant Disease",
  "Polyneuropathy In Malignant Disease",
  # 其他无法确定原发灶
  "Second Primary Malignancy",
  "Maternal Cancer In Pregnancy",
  "Cancer Cells Urine Present",
  "Malignant Neoplasm Oligoprogression",
  "Peritoneal Carcinoma Metastatic",
  # 2026-01-28 JADER完整映射
  "",
  "Carcinoma",
  "In Remission",
  "Malignant Tumour",
  "Advanced",
  "Anti-Cancer",
  "Anti-Cancer Agent",
  "Cancer Cachexia",
  "Cancer Gene Carrier",
  "Cancer Of Unknown Primary",
  "Cancer Of Unknown Primary Cause",
  "Cancer Of Unknown Primary Cause(Poorly Differentiated Adenocarcinoma",
  "Cancer(Primary Site Unknown",
  "Carcinoma Excision",
  "Carcinomatous Peritonitis Of Unknown Primary Cancer",
  "Chemotherapy",
  "Clinical",
  "Controlling Cancer Pain",
  "Decreased Blood Cells After Transplantation For Solid Cancer",
  "Egfr",
  "Epidermal Growth Factor Receptor (Egfr",
  "Gynecologic Cancer",
  "Her2 (Erbb2",
  "Her2 (Hrbb2",
  "Her2 Gene Mutation-Positive Unresectable Advanced Recurrent Non-Small Cell Lung Cancer That Has Progressed After Cancer Chemotherapy",
  "High-Dose Chemotherapy (Hdc",
  "Hypercalcaemia Due To Cancer",
  "Iatrogenic Metastasis",
  "Ld",
  "Leed Preconditioning Regimen",
  "Locally Advanced Unresectable Or Metastatic Solid Tumors",
  "Lumbar Single Metastases",
  "Management Of Cancer Pain",
  "Malignant Neoplasm Of Pelvis",
  "Metastatic Cancer",
  "Moderately Differentiated Tubular Adenocarcinoma",
  "Neoplasm Malignant (Recurrent Subclavian Lymph Node Metastasis",
  "Palliative Treatment",
  "Pathological",
  "Pelvis Metastasis",
  "Peritoneal Metastasis",
  "Poorly Differentiated Adenocarcinoma",
  "Postoperative Recurrence Of Squamous Cell Carcinoma",
  "Precordium Skin Metastasis",
  "Salvage Cbt",
  "Small?Cell Carcinoma",
  "Smarca2/Brm-Deficient Undifferentiated/Rhabdoid Carcinoma Of An Unknown Primary Site With Multiple Bone Metastases And Intraperitoneal Lymph Node Metastasis",
  "Undifferentiated Carcinoma",
  "Unknown Primary Cancer",
  "Abdominal Lymph Node Metastases",
  "Amyloidosis",
  "Ct1bn2m0",
  "Ct3n0m0",
  "Deletion",
  "Extramammary Paget Disease (L2 Bone Metastases, Lymph Node Metastases",
  "Extramammary Paget's Disease (Empd",
  "Her2-Positive Metastatic Extramammary Paget's Disease",
  "Interstitial Lung Disease",
  "Leukoencephalopathy",
  "Mcn With An Associated Invasive Carcinoma With Peritoneal Dissemination And Lymph Node Metastases",
  "Mediastinal Subcarinal Lymph Nodes",
  "Metastasis Of The Dcc",
  "Metastatic Extramammary Paget's Disease",
  "Metastatic Extramammary Paget's Disease (Empd",
  "Non-Cancer Chronic Chest Pain",
  "Nut Carcinoma",
  # 2026-01-28 JADER最终映射：原发不明癌
  "Cancer (Nos",
  "Bone Lesions Due To Solid Cancer Bone Metastasis",
  "Bone Metastasis"
)

## Hematologic（血液系统）####
## 包含：白血病、淋巴瘤、骨髓瘤等非实体瘤
organ_Hematologic <- c(
  # 2026-01-20 新增未映射PT
  "Lymphoid Leukaemia",
  "High Grade B-Cell Lymphoma",
  "Monocytic Leukaemia",
  "Non-Hodgkin's Lymphoma Transformed",
  # 浆细胞肿瘤
  "Plasma Cell Myeloma",
  "Plasma Cell Myeloma In Remission",
  "Plasma Cell Myeloma Recurrent",
  "Plasma Cell Myeloma Refractory",
  "Plasmacytoma",
  "Plasma Cell Leukaemia",
  "Plasma Cell Leukaemia In Remission",
  "Plasmablastic Lymphoma",
  "Myeloma Cast Nephropathy",
  # 慢性淋巴细胞白血病
  "Chronic Lymphocytic Leukaemia",
  "Chronic Lymphocytic Leukaemia Recurrent",
  "Chronic Lymphocytic Leukaemia Refractory",
  "Chronic Lymphocytic Leukaemia (In Remission)",
  "Chronic Lymphocytic Leukaemia Transformation",
  "B-Cell Small Lymphocytic Lymphoma",
  "B-Cell Small Lymphocytic Lymphoma Recurrent",
  "B-Cell Small Lymphocytic Lymphoma Refractory",
  "T-Cell Chronic Lymphocytic Leukaemia",
  # 慢性髓系白血病
  "Chronic Myeloid Leukaemia",
  "Chronic Myeloid Leukaemia (In Remission)",
  "Chronic Myeloid Leukaemia Recurrent",
  "Chronic Myeloid Leukaemia Transformation",
  "Philadelphia Positive Chronic Myeloid Leukaemia",
  "Blast Crisis In Myelogenous Leukaemia",
  "Myeloblastoma",
  # 骨髓增生异常/骨髓增殖性肿瘤
  "Myelodysplastic Syndrome",
  "Myeloproliferative Neoplasm",
  "Myelodysplastic Syndrome Unclassifiable",
  "Myelodysplastic Syndrome With Excess Blasts",
  "Myelodysplastic Syndrome With Multilineage Dysplasia",
  "Myelodysplastic Syndrome With Single Lineage Dysplasia",
  "Myelodysplastic Syndrome With Ringed Sideroblasts",
  "Myelodysplastic Syndrome Transformation",
  "Myeloproliferative Disorder",
  # 急性淋巴细胞白血病
  "Acute Lymphocytic Leukaemia",
  "Acute Lymphocytic Leukaemia Recurrent",
  "Acute Lymphocytic Leukaemia Refractory",
  "Acute Lymphocytic Leukaemia (In Remission)",
  "Philadelphia Positive Acute Lymphocytic Leukaemia",
  "B Precursor Type Acute Leukaemia",
  "B-Cell Type Acute Leukaemia",
  "T-Cell Type Acute Leukaemia",
  "Precursor T-Lymphoblastic Leukaemia Acute",
  "Myeloblastoma",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia Refractory",
  "Precursor T-Lymphoblastic Lymphoma/Leukaemia Recurrent",
  "Precursor B-Lymphoblastic Lymphoma",
  "Precursor B-Lymphoblastic Lymphoma Recurrent",
  "Precursor B-Lymphoblastic Lymphoma Refractory",
  "Acute Biphenotypic Leukaemia",
  "Acute Bilineal Leukaemia",
  "Acute Undifferentiated Leukaemia",
  "Mature B-Cell Type Acute Leukaemia",
  # B细胞淋巴瘤
  "B-Cell Lymphoma",
  "B-Cell Lymphoma Refractory",
  "B-Cell Lymphoma Recurrent",
  "B-Cell Unclassifiable Lymphoma High Grade",
  "B-Cell Unclassifiable Lymphoma Low Grade",
  "B-Cell Lymphoma Unclassifiable",
  "High-Grade B-Cell Lymphoma",
  "High Grade B-Cell Lymphoma Refractory",
  "High Grade B-Cell Lymphoma Recurrent",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma Recurrent",
  "High Grade B-Cell Lymphoma Burkitt-Like Lymphoma Refractory",
  "Primary Mediastinal Large B-Cell Lymphoma",
  "Primary Mediastinal Large B-Cell Lymphoma Refractory",
  "Primary Mediastinal Large B-Cell Lymphoma Recurrent",
  "Primary Gastrointestinal Follicular Lymphoma",
  # 套细胞淋巴瘤
  "Mantle Cell Lymphoma",
  "Mantle Cell Lymphoma Recurrent",
  "Mantle Cell Lymphoma Refractory",
  # 非霍奇金淋巴瘤
  "Non-Hodgkin's Lymphoma",
  "Non-Hodgkin's Lymphoma Recurrent",
  "Non-Hodgkin's Lymphoma Refractory",
  "Non-Hodgkin's Lymphoma Unspecified Histology Indolent",
  "Non-Hodgkin's Lymphoma Unspecified Histology Aggressive",
  "Non-Hodgkin's Lymphoma Unspecified Histology Aggressive Refractory",
  "Non-Hodgkin's Lymphoma Metastatic",
  "Non-Hodgkin's Lymphoma Transformed Recurrent",
  # 霍奇金淋巴瘤
  "Hodgkin's Disease",
  "Hodgkin's Disease Nodular Sclerosis",
  "Hodgkin's Disease Refractory",
  "Hodgkin's Disease Recurrent",
  "Hodgkin's Disease Mixed Cellularity",
  "Hodgkin's Disease Mixed Cellularity Recurrent",
  "Hodgkin's Disease Mixed Cellularity Refractory",
  "Hodgkin's Disease Lymphocyte Predominance Type",
  "Hodgkin's Disease Lymphocyte Predominance",
  "Hodgkin's Disease Lymphocyte Depletion Type",
  "Hodgkin's Disease Lymphocyte Depletion",
  "Hodgkin's Disease Lymphocyte Depletion Type Refractory",
  "Hodgkin's Disease Nodular Sclerosis Recurrent",
  "Hodgkin's Disease Unclassifiable",
  "Nodular Lymphocyte Predominant Hodgkin Lymphoma",
  # 滤泡性淋巴瘤
  "Follicular Lymphoma",
  "Follicular Lymphoma Refractory",
  "Follicular Lymphoma Recurrent",
  "Follicle Centre Lymphoma, Follicular",
  "Follicle Centre Lymphoma Diffuse Small Cell Lymphoma",
  "Follicle Centre Lymphoma Diffuse Small Cell Lymphoma Recurrent",
  # T细胞淋巴瘤
  "T-Cell Lymphoma",
  "T-Cell Lymphoma Refractory",
  "T-Cell Lymphoma Recurrent",
  "Cutaneous T-Cell Lymphoma",
  "Cutaneous T-Cell Lymphoma Refractory",
  "Cutaneous T-Cell Lymphoma Recurrent",
  "Peripheral T-Cell Lymphoma Unspecified",
  "Peripheral T-Cell Lymphoma Unspecified Refractory",
  "Peripheral T-Cell Lymphoma Unspecified Recurrent",
  "Angioimmunoblastic T-Cell Lymphoma",
  "Angioimmunoblastic T-Cell Lymphoma Refractory",
  "Angioimmunoblastic T-Cell Lymphoma Recurrent",
  "Adult T-Cell Lymphoma/Leukaemia",
  "Adult T-Cell Lymphoma/Leukaemia Recurrent",
  "Adult T-Cell Lymphoma/Leukaemia Refractory",
  "Enteropathy-Associated T-Cell Lymphoma",
  "Intestinal T-Cell Lymphoma",
  "Intestinal T-Cell Lymphoma Recurrent",
  "Intestinal T-Cell Lymphoma Refractory",
  "Hepatosplenic T-Cell Lymphoma",
  "Angiocentric Lymphoma",
  "Angiocentric Lymphoma Recurrent",
  "Natural Killer-Cell Lymphoblastic Lymphoma",
  "Extranodal Nk/T-Cell Lymphoma, Nasal Type",
  # 间变性大细胞淋巴瘤
  "Anaplastic Large-Cell Lymphoma",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types Refractory",
  "Anaplastic Large Cell Lymphoma T- And Null-Cell Types Recurrent",
  "Anaplastic Lymphoma Kinase Gene Mutation",
  "Anaplastic Lymphoma Kinase Gene And Nucleophosmin Gene Fusion Overexpression",
  # 边缘区淋巴瘤
  "Marginal Zone Lymphoma",
  "Marginal Zone Lymphoma Recurrent",
  "Marginal Zone Lymphoma Refractory",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type)",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type) Recurrent",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type) Refractory",
  "Extranodal Marginal Zone B-Cell Lymphoma (Balt Type)",
  "Nodal Marginal Zone B-Cell Lymphoma",
  "Nodal Marginal Zone B-Cell Lymphoma Recurrent",
  "Nodal Marginal Zone B-Cell Lymphoma Refractory",
  "Splenic Marginal Zone Lymphoma",
  "Splenic Marginal Zone Lymphoma Recurrent",
  "Splenic Marginal Zone Lymphoma Refractory",
  # 伯基特淋巴瘤
  "Burkitt's Lymphoma",
  "Burkitt's Lymphoma Recurrent",
  "Burkitt's Lymphoma Refractory",
  "Burkitt's Leukaemia",
  # 华氏巨球蛋白血症/淋巴浆细胞淋巴瘤
  "Lymphoplasmacytoid Lymphoma/Immunocytoma",
  "Lymphoplasmacytoid Lymphoma/Immunocytoma Recurrent",
  "Lymphoplasmacytoid Lymphoma/Immunocytoma Refractory",
  # 其他淋巴瘤
  "Lymphoma",
  "Lymphoma Transformation",
  "Lymphoma Aids Related",
  "Gastrointestinal Lymphoma",
  "Ocular Lymphoma",
  "Lymphocytic Lymphoma",
  "Disseminated Large Cell Lymphoma",
  "Primary Effusion Lymphoma",
  "Primary Cardiac Lymphoma",
  "Double Hit Lymphoma",
  "Grey Zone Lymphoma",
  "Triple Hit Lymphoma",
  "Composite Lymphoma",
  "Immunoblastic Lymphoma",
  "Lymphoma Cutis",
  "Metastatic Lymphoma",
  "Leukaemic Lymphoma",
  "Epstein-Barr Virus Associated Lymphoma",
  "Malignant Lymphoma Unclassifiable High Grade",
  # 蕈样肉芽肿
  "Mycosis Fungoides",
  "Mycosis Fungoides Refractory",
  # 其他白血病
  "Leukaemia",
  "Leukaemia Recurrent",
  "Leukaemia In Remission",
  "Lymphocytic Leukaemia",
  "Myeloid Leukaemia",
  "Myeloid Leukaemia In Remission",
  "Acute Leukaemia",
  "Acute Leukaemia In Remission",
  "Chronic Leukaemia",
  "Chronic Leukaemia In Remission",
  "Chronic Myelomonocytic Leukaemia",
  "Chronic Myelomonocytic Leukaemia (In Remission)",
  "Juvenile Chronic Myelomonocytic Leukaemia",
  "Hairy Cell Leukaemia",
  "Hairy Cell Leukaemia Recurrent",
  "T-Cell Prolymphocytic Leukaemia",
  "B-Cell Prolymphocytic Leukaemia",
  "Prolymphocytic Leukaemia",
  "Natural Killer-Cell Leukaemia",
  "Chronic Eosinophilic Leukaemia",
  "Eosinophilic Leukaemia",
  "Mastocytic Leukaemia",
  "Leukaemia Cutis",
  "Leukaemia Basophilic",
  "Leukaemia Granulocytic",
  "Leukaemia Monocytic",
  "Monocytic Leukaemia In Remission",
  "Lymphoid Leukaemia (In Remission)",
  "Leukaemia Recurrent Prophylaxis",
  "Lineage Switch Leukaemia",
  "Aleukaemic Leukaemia",
  # 血液系统恶性肿瘤（泛称）
  "Haematological Malignancy",
  # 2026-01-28 JADER新增（日英名称不同的PT）
  "Multiple Myeloma",
  "Relapsed/Refractory Multiple Myeloma",
  "Relapsed Or Refractory Multiple Myeloma",
  "Refractory Multiple Myeloma",
  "Multiple Myeloma Recurrent",
  "Chronic Myeloid Leukemia",
  "Chronic Myelogenous Leukemia",
  "Adult T-Cell Leukemia-Lymphoma",
  "Adult T-Cell Leukemia/Lymphoma",
  "Adult T-Cell Leukemia Lymphoma",
  "Adult T-Cell Leukemia",
  "Acute Lymphoblastic Leukemia",
  "Acute Lymphocytic Leukemia",
  "Malignant Lymphoma",
  "Hodgkin's Lymphoma",
  "Classical Hodgkin Lymphoma",
  "Lymphoplasmacytic Lymphoma",
  "Primary Central Nervous System Lymphoma",
  # 2026-01-28 JADER新增（第二批）
  "Malt Type",
  "Cd22 Positive Acute Lymphocytic Leukaemia",
  "Cd22-Positive Acute Lymphocytic Leukaemia",
  "5q Minus Myelodysplastic Syndrome",
  "Acute Lymphoblastic Leukaemia",
  "T-Cell Acute Lymphocytic Leukaemia",
  "Follicular B-Cell Non-Hodgkin's Lymphoma",
  "T-Lymphoblastic Lymphoma",
  # 2026-01-28 JADER完整映射（血液系统）
  "Relapsed Or Refractory Multiple Myeloma (When Standard Therapy Is Difficult",
  "Relapsed Or Refractory Multiple Myeloma (Only When Standard Treatment Is Difficult",
  "Relapsed And Refractory Multiple Myeloma (Rrmm",
  "Relapsed Or Refractory Multiple Myeloma (Limited To Cases Where The Standard Of Care Is Difficult",
  "Treatment Of Multiple Myeloma",
  "Multiple Myeloma (Mm",
  "Igg-Λ Relapsed Refractory Multiple Myeloma(Rrmm",
  "Numbness Due To Multiple Myeloma (Mm",
  "Multiple Myeloma (Refractory",
  "Recurrent Multiple Myeloma / Refractory Multiple Myeloma",
  "Relapsed And Refractory Multiple Myeloma (Limited To Cases Where The Standard Of Care Is Difficult",
  "Relapsed Or Refractory Multiple Myeloma (Only If Standard Treatment Was Difficult",
  "Relapsed Refractory Multiple Myeloma",
  "Relapsed/Refractory Multiple Myeloma (Only When Standard Treatment Is Difficult",
  "Relapsed/Refractory Multiple Myeloma 　",
  "Relapsed/Refractory Multiple Myeloma（Rrmm",
  "Symptomatic Myeloma, Igg Type",
  "Treatment For Multiple Myeloma",
  "Treatment Of Multiple Myeloma/Relapsed And Refractory Multiple Myeloma",
  "Cd22-Possitive Acute Lymphocytic Leukaemia",
  "Chronic Phase Chronic Myeloid Leukaemia",
  "Extranodal Marginal Zone B-Cell Lymphoma (Malt Type",
  "Folliculotropic Mycosis Fungoides",
  "Lymphoblastic Leukemia",
  "Nodular (Follicular",
  "Acute Lymphocytic Leukaemia (All",
  "All (Acute Lymphocytic Leukaemia",
  "Acute Lymphocytic Leukaemia(B-All",
  "Acute Lymphocytic Leukaemia Recurrence",
  "Acute Lymphocytic Leukaemia/Philadelphia Negative Acute Lymphocytic Leukaemia",
  "Acute Megakaryoblastic Leukemia (Amkl",
  "Acute Adult T-Cell Leukemia",
  "Acute Leukemia",
  "Acute Type Of Adult T-Cell Leukemia/Lymphoma",
  "Acute Unclassified Leukemia",
  "Adult T-Cell Leukemia-Lymphoma, Acute Type",
  "Adult T-Cell Leukemia/Lymphoma (Atll",
  "Adult T-Cell Lymphoma",
  "Anaplastic Large Cell Lymphoma",
  "Angioimmunoblastic T Cell Lymphoma",
  "Angioimmunoblastic T-Cell Lymphoma (Aitl",
  "Angioimmunoblastic T-Cell Lymphoma Complicated By Haemophagocytic Syndrome",
  "Burkitt Lymphoma",
  "Burkitt's Malignant Lymphoma",
  "Cd22 Positive Acute Lymphocytic Leukaemia Refractory And",
  "Cd22 Positive Acute Lymphocytic Leukaemia Refractory/Transplant Ineligible Acute Lymphocytic Leukaemia (All",
  "Cd22-Positive Acute Lymphocytic Leukaemia (Recurrence",
  "Cd22-Positive Acute Lymphocytic Leukemia (Intractable",
  "Cd30-Positive Primary Effusion Lymphoma-Like Adult T-Cell Leukemia/Lymphoma",
  "Central Nerve Invasion Of T-Lymphoblastic Lymphoma",
  "Central Recurrence Of Intra-Ocular Lymphoma/ Brain Malignant Lymphoma",
  "Childhood Acute Lymphocytic Leukaemia",
  "Chronic Active Epstein-Bar Virus Infection",
  "Chronic Lymphocytic Leukemia",
  "Chronic Myeloid Leukaemia (Cml",
  "Chronic Myeloid Leukaemia (Cml, Initial Episode",
  "Chronic Myeloid Leukaemia Chronic Phase",
  "Chronic Myeloid Leukaemia Reoccurred",
  "Chronic Myeloid Leukaemia(Cml",
  "Chronic Myeloid Leukaemia, First Time Therapy",
  "Chronic Myeloid Leukaemia/Cml-Ap First-Line",
  "Chronic Myeloid Leukemia Blast Phase ( Lymphoid Crisis",
  "Chronic Myeloid Leukemia Transformation",
  "Chronic Neutrophilic Leukemia",
  "Chronic Phase Chronic Myeloid Leukaemia (Cml-Cp",
  "Classic Hodgkin Lymphoma",
  "Classical Hodgkin Lymphoma, Nodular Sclerosis,",
  "Cml(Chronic Myeloid Leukemia",
  "Consolidation For Acute Lymphocytic Leukaemia",
  "Diffuse Large B Cell Lymphoma",
  "Diffuse Large B Cell Lymphoma (Dlbcl",
  "Follicular Lymphoma (Fl",
  "Follicular Lymphoma(Fl",
  "Follicular Lymphoma（Fl）",
  "Gastric Burkitt Lymphoma",
  "Hairy Cell Leukemia",
  "Hematological Malignancy",
  "Hodgkin's Lymphoma, Mixed Cellularity Type,",
  "Indolent B-Cell Lymphoma＋Myelofibrosis",
  "Induction Of Remission Of Leukaemia",
  "Intraocular Malignant Lymphoma",
  "Intravascular Large B-Cell Lymphoma",
  "Intravascular Large B-Cell Lymphoma (Ivlbcl",
  "Intravascular Large Cell Lymphoma",
  "Intravascular Lymphoma",
  "Isochromosome 1qを伴うEbv陽性Burkitt Leukemia/Lymphoma",
  "Large B-Cell Lymphoma",
  "Leukaemia With Central Nervous System Infiltration",
  "Leukemia Remission Induction",
  "Lymphoblastic Blast Crisis Of Chronic Myelocytic Leukemia",
  "Lymphoplasmacytic Lymphoma/Waldenstorm's Macroglobulinemia(Lpl/Wm",
  "Lymphoplasmacytic Type Lymphoma",
  "Malignant Lymphoma (Dlbcl",
  "Mature B-Cell Acute Lymphocytic Leukemia",
  "Mds Overt Leukemia",
  "Metastatic Mantle Cell Lymphoma",
  "Mixed-Cellularity Classical Hodgkin Lymphoma, Cs Iia, And Ps 0",
  "Mixed-Phenotype Acute Leukaemia",
  "Mixed-Phenotype Acute Leukemia (Mpal",
  "Myelodysplasla-Related Changes",
  "Myelodysplastic Syndrome (Mds",
  "Myelodysplastic Syndrome-Mds",
  "Myelodysplastic/Myeloproliferative Disease",
  "Myelofibrosis Associated With Myelodysplastic Syndrome",
  "Myelofibrosis Progressed From Chronic Myelomonocytic Leukaemia",
  "Myeloid Leukemia Associated With Down Syndrome",
  "Nodal Epstein-Barr Virus-Positive T/Nk-Cell Lymphoma",
  "Non-Hodgkin Lymphoma、Diffuse Small B Cell Lymphoma",
  "Non‐Hodgkin Lymphoma",
  "Pancreatic Malt Lymphoma",
  "Pediatric Anaplastic Large-Cell Lymphoma",
  "Peripheral T-Cell Lymphoma",
  "Ph Positive Acute Lymphocytic Leukaemia",
  "Ph-Like Acute Lymphoblastic Leukaemia",
  "Ph-Positive Acute Lymphocytic Leukaemia",
  "Philadelphia Chromosome Negative B-Cell Acute Lymphoblastic Leukaemia",
  "Philadelphia Chromosome Negative B-Cell Acute Lymphocytic Leukemia",
  "Philadelphia Chromosome-Positive Acute Lymphocytic Leukaemia",
  "Philadelphia Negative, B-Cell Type Acute Lymphoblastic Leukaemia",
  "Philadelphia Positive Acute Lymphocytic Leukaemia (A Relapse",
  "Plasma Cell Leukemia",
  "Pleural Efusion With Lymphoma Cell Infltration",
  "Post-Unclassifiable Myeloproliferative Disorder Myelofibrosis",
  "Precursor B Acute Lymphocytic Leukaemia",
  "Precursor B-Cell Acute Lymphoblastic Leukaemia",
  "Primary Central Nervous Symptom Lymphoma",
  "Primary Central Nervous System Lymphoma (Pcnsl",
  "Primary Central Malignant Lymphoma",
  "Primary Central Nervous System Low-Grade B-Cell Lymphoma",
  "Primary Central Nervous System Malignant Lymphoma",
  "Primary Central Nervous System Malignant Lymphoma (Pcnsl",
  "Primary Chronic Myeloid Leukaemia",
  "Primary Cml (Chronic Myelogenous Leukemia）",
  "Primary Leptomeningeal B-Cell Lymphoma",
  "Primary Malignant Lymphoma Of The Central Nervous System",
  "Primary Malignant Lymphoma Of The Central Nervous System (Pcnsl",
  "Recurrence Of Philadelphia Chromosome Negative B-Cell Acute Lymphocytic Leukaemia",
  "Recurrent Acute Lymphocytic Leukaemia",
  "Recurrent Cd22-Positive Acute Lymphocytic Leukaemia",
  "Recurrent Or Refractory B-Cell Acute Lymphocytic Leukaemia (B Precursor Type Acute Leukaemia",
  "Recurrent Refractory Acute Lymphocytic Leukaemia(R/Rall",
  "Recurrent/Refractory Acute Lymphocytic Leukaemia",
  "Refractory Acute Lymphoblastic Leukemia",
  "Refractory Cd22-Positive B-Cell Acute Lymphoblastic Leukemia",
  "Refractory, Cd22 Positive Acute Lymphocytic Leukaemia",
  "Relapsed Acute Lymphoblastic Leukemia",
  "Relapsed Follicular Lymphoma",
  "Relapsed Or Refractory Adult T-Cell Leukemia-Lymphoma",
  "Relapsed Refractory Acute Lymphoblastic Leukemia",
  "Relapsed/Refractory Adult T-Cell Leukemia-Lymphoma",
  "Secondary Itp With Mantle Cell Lymphoma",
  "Secondary Organising Pneumonia Of Myelodysplastic Syndrome (Mds",
  "T Cell Large Granular Lymphocyte (T-Lgl",
  "T Lymphoblastic Lymphoma (Murphy",
  "T-Cell Acute Lymphoblastic Leukaemia",
  "T-Cell Lymphoblastic Leukemia/ Lymphoma",
  "T-Cell Lymphoblastic Lymphoma",
  "T-Cell Lymphoblastic Lymphoma (T-All",
  "T-Cell/Histiocyte Rich B-Cell Lymphoma",
  "Treatment For Acute Lymphoblastic Leukaemia",
  "Treatment For Chronic Myeloid Leukaemia",
  "Waldenstrom's Macroglobulinaemia",
  "Waldenstrom's Macroglobulinemia/Lymphoplasmacytic Lymphoma (Wm/Lpl",
  "Acml",
  "Advanced Extranodal Nk/T-Cell Lymphoma, Nasal Type",
  "Extranodal Nk/T-Cell Lymphoma, Nasal Type (Enkl",
  "Mf",
  # 2026-01-28 JADER最终映射
  "Asymptomatic Chronic Lymphocytic Leukemia With Progressive Splenomegaly",
  "Atl (Adult T-Cell Leukemia",
  "B Cell Precursor Acute Lymphoblastic Leukemia",
  "B-Cell Acute Lymphoblastic Leukaemia",
  "B-Cell Acute Lymphoblastic Leukemia",
  "B-Cell Acute Lymphocytic Leukaemia",
  "B-Cell Lymphoblastic Lymphoma (",
  "B-Cell Non Hodgkin Lymphoma",
  "B-Cell Non-Hodgkin's Lymphoma",
  "B-Cell Precursor Acute Lymphoblastic Leukemia",
  "B-Chronic Lymphocytic Leukemia",
  "B-Precursor Cell Acute Lymphoblastic Leukaemia",
  "Bcma Against Multiple Myeloma",
  "Bone Lesions Due To Multiple Myeloma",
  "Lymphoplasmacytic Lymphoma/Waldenstorm's Macroglobulinemia(Lpl/Wm",
  # 2026-01-28 JADER最终映射：华氏巨球蛋白血症变体（不同撇号/空格格式）
  "Lymphoplasmacytic Lymphoma/Waldenstorm's Macroglobulinemia(Lpl/Wm",
  "Waldenstrom's Macroglobulinemia/Lymphoplasmacytic Lymphoma (Wm/Lpl"
)

## Rare（稀有肿瘤）####
## 包含：能确认原发灶但不属于前面大类的稀有肿瘤
## 2026-01-20：合并了Ampulla of Vater, Penis/Male Genital, Peritoneum, Vulva/Vagina
organ_Rare <- c(
  # 肾上腺
  "Adrenal Gland Cancer",
  "Adrenal Gland Cancer Metastatic",
  # 胎盘
  "Malignant Neoplasm Of Placenta",
  # 脾脏
  "Splenic Neoplasm Malignancy Unspecified",
  # 心脏
  "Cardiac Neoplasm Malignant",
  # 胸部
  "Malignant Neoplasm Of Thorax",
  # 腹壁
  "Abdominal Wall Neoplasm Malignant",
  # 颌面部（成釉细胞）
  "Ameloblastoma",
  "Ameloblastic Carcinoma",
  # 其他稀有（非肿瘤相关项已移除）
  "Anti-Melanoma Differentiation-Associated Protein 5 Antibody Positive",
  "Kaposi's Varicelliform Eruption",
  # 睾丸/生殖细胞肿瘤（原Testis类别并入）
  "Germ Cell Cancer",
  "Germ Cell Cancer Metastatic",
  "Extragonadal Primary Embryonal Carcinoma",
  "Testicular Cancer Metastatic",
  # 壶腹部（原Ampulla_of_Vater类别并入）
  "Malignant Neoplasm Of Ampulla Of Vater",
  "Malignant Neoplasm Papilla Of Vater",
  # 腹膜（原Peritoneum类别并入）
  "Metastases To Peritoneum",
  "Peritoneal Cancer Index",
  "Peritoneal Carcinoma",
  # 外阴/阴道（原Vulva_Vagina类别并入）
  "Vulval Cancer",
  "Vulval Cancer Recurrent",
  "Vulval Cancer Metastatic",
  "Squamous Cell Carcinoma Of The Vulva",
  "Vulvar Adenocarcinoma",
  "Vulvar Basal Cell Carcinoma",
  "Vaginal Cancer",
  "Vaginal Cancer Recurrent",
  "Vaginal Adenocarcinoma",
  "Squamous Cell Carcinoma Of The Vagina",
  "Genital Neoplasm Malignant Female",
  "Female Reproductive Tract Carcinoma In Situ",
  "Vaginal Cancer Metastatic",
  # 阴茎/男性生殖器（原Penis_Male_Genital类别并入）
  "Penile Cancer",
  "Penile Squamous Cell Carcinoma",
  "Penis Carcinoma",
  "Penis Carcinoma Metastatic",
  "Carcinoma In Situ Of Penis",
  "Genital Cancer Male",
  "Scrotal Cancer",
  "Epididymal Cancer",
  "Malignant Neoplasm Of Spermatic Cord"
)

## Eye（眼）####
organ_Eye <- c(
  "Retinoblastoma",
  "Malignant Neoplasm Of Choroid",
  "Malignant Neoplasm Of Eye",
  "Ocular Lymphoma",
  "Malignant Neoplasm Of Conjunctiva",
  "Malignant Neoplasm Of Lacrimal Duct",
  "Malignant Neoplasm Of Lacrimal Gland",
  "Malignant Neoplasm Of Retina",
  "Malignant Neoplasm Of Orbit",
  "Malignant Neoplasm Of Eyelid",
  "Malignant Neoplasm Of Cornea",
  "Congenital Retinoblastoma",
  "Extraocular Retinoblastoma",
  "Carcinoma In Situ Of Eye",
  "Ocular Cancer Metastatic",
  # 移除Gastrointestinal Melanoma和Genitourinary Melanoma（不属于眼）
  "Gastrointestinal Melanoma",
  "Genitourinary Melanoma",
  # 2026-01-20 新增未映射PT
  "Ocular Cancer"
)


# 映射函数定义####

## 主映射函数：同时支持 assigned_type 和 first_pt####
## 参数说明：
##   - first_pt: PT名称（首字母大写格式）
##   - assigned_type: TCGA癌种代码（可选，默认为NA）
## 匹配优先级：
##   1. 如果 assigned_type 有值且非 "Others"，优先使用 tcga_to_organ 映射
##   2. 否则使用 first_pt 在 organ_xxx 向量中匹配
## 分类参考：泛癌种分类层级图
map_to_primary_organ <- function(first_pt, assigned_type = NA_character_) {
  
  # 优先级1：通过 assigned_type（TCGA代码）映射
  if (!is.na(assigned_type) && assigned_type != "Others" && assigned_type != "") {
    organ <- tcga_to_organ[[assigned_type]]
    if (!is.null(organ)) {
      return(organ)
    }
  }
  
  # 优先级2：通过 first_pt 在 organ_xxx 向量中匹配
  if (!is.na(first_pt) && first_pt != "") {
    # 规范化PT格式为Title Case（与organ_xxx向量保持一致）
    first_pt_normalized <- stringr::str_to_title(first_pt)
    
    if (first_pt_normalized %in% organ_Lung) return("Lung")
    if (first_pt_normalized %in% organ_Pleura) return("Pleura")
    if (first_pt_normalized %in% organ_Kidney) return("Kidney")
    if (first_pt_normalized %in% organ_Bladder) return("Bladder")
    if (first_pt_normalized %in% organ_Prostate) return("Prostate")
    if (first_pt_normalized %in% organ_Cervix) return("Cervix")
    if (first_pt_normalized %in% organ_Uterus) return("Uterus")
    if (first_pt_normalized %in% organ_Breast) return("Breast")
    if (first_pt_normalized %in% organ_Ovary) return("Ovary")
    if (first_pt_normalized %in% organ_Bowel) return("Bowel")
    if (first_pt_normalized %in% organ_Esophagus_stomach) return("Esophagus/stomach")
    if (first_pt_normalized %in% organ_Liver) return("Liver")
    if (first_pt_normalized %in% organ_Biliary_tract) return("Biliary tract")
    if (first_pt_normalized %in% organ_Pancreas) return("Pancreas")
    if (first_pt_normalized %in% organ_Head_neck) return("Head/neck")
    if (first_pt_normalized %in% organ_Skin) return("Skin")
    if (first_pt_normalized %in% organ_Brain) return("Brain")
    if (first_pt_normalized %in% organ_PNS) return("PNS")
    if (first_pt_normalized %in% organ_Bone) return("Bone")
    if (first_pt_normalized %in% organ_Soft_tissue) return("Soft tissue")
    if (first_pt_normalized %in% organ_Thymus) return("Thymus")
    if (first_pt_normalized %in% organ_Thyroid) return("Thyroid")
    if (first_pt_normalized %in% organ_CUP) return("CUP")
    if (first_pt_normalized %in% organ_Hematologic) return("Hematologic")
    if (first_pt_normalized %in% organ_Rare) return("Rare")
    if (first_pt_normalized %in% organ_Eye) return("Eye")
  }
  
  # 未匹配到的返回NA
  return(NA_character_)
}

## 仅通过 assigned_type 映射（快速版本）####
map_tcga_to_organ <- function(assigned_type) {
  if (is.na(assigned_type) || assigned_type == "Others" || assigned_type == "") {
    return(NA_character_)
  }
  organ <- tcga_to_organ[[assigned_type]]
  if (is.null(organ)) {
    return(NA_character_)
  }
  return(organ)
}

## 仅通过 first_pt 映射（原始版本，向后兼容）####
map_pt_to_organ <- function(first_pt) {
  if (is.na(first_pt) || first_pt == "") {
    return(NA_character_)
  }
  if (first_pt %in% organ_Lung) return("Lung")
  if (first_pt %in% organ_Pleura) return("Pleura")
  if (first_pt %in% organ_Kidney) return("Kidney")
  if (first_pt %in% organ_Bladder) return("Bladder")
  if (first_pt %in% organ_Prostate) return("Prostate")
  if (first_pt %in% organ_Cervix) return("Cervix")
  if (first_pt %in% organ_Uterus) return("Uterus")
  if (first_pt %in% organ_Breast) return("Breast")
  if (first_pt %in% organ_Ovary) return("Ovary")
  if (first_pt %in% organ_Bowel) return("Bowel")
  if (first_pt %in% organ_Esophagus_stomach) return("Esophagus/stomach")
  if (first_pt %in% organ_Liver) return("Liver")
  if (first_pt %in% organ_Biliary_tract) return("Biliary tract")
  if (first_pt %in% organ_Pancreas) return("Pancreas")
  if (first_pt %in% organ_Head_neck) return("Head/neck")
  if (first_pt %in% organ_Skin) return("Skin")
  if (first_pt %in% organ_Brain) return("Brain")
  if (first_pt %in% organ_PNS) return("PNS")
  if (first_pt %in% organ_Bone) return("Bone")
  if (first_pt %in% organ_Soft_tissue) return("Soft tissue")
  if (first_pt %in% organ_Thymus) return("Thymus")
  if (first_pt %in% organ_Thyroid) return("Thyroid")
  if (first_pt %in% organ_CUP) return("CUP")
  if (first_pt %in% organ_Hematologic) return("Hematologic")
  if (first_pt %in% organ_Rare) return("Rare")
  if (first_pt %in% organ_Eye) return("Eye")
  return(NA_character_)
}

## 获取所有已定义的器官列表####
get_organ_list <- function() {
  c(
    "Lung", "Pleura", "Kidney", "Bladder", "Prostate", "Cervix",
    "Uterus", "Breast", "Ovary", "Bowel",
    "Esophagus/stomach", "Liver", "Biliary tract",
    "Pancreas", "Head/neck", "Skin", "Brain", "PNS", "Bone",
    "Soft tissue", "Thymus", "Thyroid", "CUP", "Hematologic",
    "Rare", "Eye"
  )
}

## 获取映射规则汇总表####
get_mapping_summary <- function() {
  organ_list <- get_organ_list()
  organ_vars <- c(
    "organ_Lung", "organ_Pleura", "organ_Kidney", "organ_Bladder",
    "organ_Prostate", "organ_Cervix", "organ_Uterus", "organ_Breast",
    "organ_Ovary", "organ_Bowel",
    "organ_Esophagus_stomach", "organ_Liver", "organ_Biliary_tract",
    "organ_Pancreas", "organ_Head_neck", "organ_Skin", "organ_Brain",
    "organ_PNS", "organ_Bone", "organ_Soft_tissue", "organ_Thymus",
    "organ_Thyroid", "organ_CUP", "organ_Hematologic",
    "organ_Rare", "organ_Eye"
  )
  
  pt_counts <- sapply(organ_vars, function(v) length(get(v)))
  
  data.frame(
    organ = organ_list,
    n_pts = pt_counts,
    stringsAsFactors = FALSE
  )
}

# 第二层映射函数：Histologic Category（组织学分类）####

## 主映射函数：同时支持 PT、assigned_type 和 primary_organ####
## 参数说明：
##   - first_pt: PT名称（首字母大写格式）
##   - assigned_type: TCGA癌种代码（可选，默认为NA）
##   - primary_organ: 第一层器官分类结果（可选，默认为NA）
## 匹配优先级：
##   1. 首先使用 first_pt 在组织学分类PT向量中匹配
##   2. 如果 assigned_type 有值，使用 tcga_to_histology 映射
##   3. 最后使用 primary_organ 在 organ_to_histology 中获取默认分类
map_to_histology <- function(first_pt, assigned_type = NA_character_, primary_organ = NA_character_) {
  
  # 优先级1：通过 first_pt 在组织学分类PT向量中匹配
  if (!is.na(first_pt) && first_pt != "") {
    # 明确的组织学类型（优先匹配）
    if (first_pt %in% histology_Adenocarcinoma) return("Adenocarcinoma")
    if (first_pt %in% histology_Squamous_Cell_Carcinoma) return("Squamous Cell Carcinoma")
    if (first_pt %in% histology_Urothelial_Carcinoma) return("Urothelial Carcinoma")
    if (first_pt %in% histology_Hepatocellular_Carcinoma) return("Hepatocellular Carcinoma")
    if (first_pt %in% histology_Renal_Cell_Carcinoma) return("Renal Cell Carcinoma")
    if (first_pt %in% histology_Thyroid_Carcinoma) return("Thyroid Carcinoma")
    if (first_pt %in% histology_Basal_Cell_Carcinoma) return("Basal Cell Carcinoma")
    if (first_pt %in% histology_Small_Cell_Carcinoma) return("Small Cell Carcinoma")
    if (first_pt %in% histology_Melanoma) return("Melanoma")
    if (first_pt %in% histology_Sarcoma) return("Sarcoma")
    if (first_pt %in% histology_Glioma) return("Glioma")
    if (first_pt %in% histology_Neuroblastoma) return("Neuroblastoma")
    if (first_pt %in% histology_Lymphoma) return("Lymphoma")
    if (first_pt %in% histology_Leukemia) return("Leukemia")
    if (first_pt %in% histology_Plasma_Cell_Neoplasm) return("Plasma Cell Neoplasm")
    if (first_pt %in% histology_MDS_MPN) return("MDS/MPN")
    if (first_pt %in% histology_Germ_Cell_Tumor) return("Germ Cell Tumor")
    if (first_pt %in% histology_Neuroendocrine_Tumor) return("Neuroendocrine Tumor")
    if (first_pt %in% histology_Embryonal_Tumor) return("Embryonal Tumor")
    if (first_pt %in% histology_Mesothelioma) return("Mesothelioma")
    if (first_pt %in% histology_Thymoma) return("Thymoma")
    if (first_pt %in% histology_Other_Carcinoma) return("Other Carcinoma")
    # Unknown放最后（无法从PT确定组织学类型）
    if (first_pt %in% histology_Unknown) return("Unknown")
  }
  
  # 优先级2：通过 assigned_type（TCGA代码）映射
  if (!is.na(assigned_type) && assigned_type != "Others" && assigned_type != "") {
    histology <- tcga_to_histology[[assigned_type]]
    if (!is.null(histology)) {
      return(histology)
    }
  }
  
  # 注意：已删除原优先级3（organ_to_histology映射）
  # 因为所有PT+TCGA代码的组合已完整覆盖，不需要额外兜底
  
  # 未匹配到的返回NA
  return(NA_character_)
}

# PT到Histologic_category映射（严格按照用户映射表Results/癌种映射.xlsx）####
# 说明：将PT直接映射到用户定义的Histologic_category
# 这个函数用于没有TCGA代码的情况
# 格式：histcat_xxx <- c("PT1", "PT2", ...)

## Bowel器官的Histologic_category####
histcat_Small_bowel_cancer <- c(
  "Small Intestine Carcinoma",
  "Small Intestine Adenocarcinoma",
  "Small Intestine Carcinoma Metastatic",
  "Small Intestine Carcinoma Recurrent"
)

histcat_Colorectal_cancer <- c(
  "Colorectal Cancer",
  "Colorectal Cancer Metastatic",
  "Colorectal Cancer Recurrent",
  "Colorectal Adenocarcinoma",
  "Colorectal Carcinoma",
  "Colon Cancer",
  "Colon Cancer Metastatic",
  "Colon Cancer Recurrent",
  "Rectal Cancer",
  "Rectal Cancer Metastatic",
  "Rectal Cancer Recurrent",
  "Rectosigmoid Cancer",
  "Rectosigmoid Cancer Recurrent",
  "Rectosigmoid Cancer Metastatic",
  "Hereditary Non-Polyposis Colorectal Cancer Syndrome",
  "Her2 Positive Colorectal Cancer",
  "Recurrent N-Ras Mutation-Positive Colorectal Carcinoma"
)

histcat_Appendiceal_cancer <- c(
  "Appendix Cancer",
  "Appendix Cancer Metastatic",
  "Adenocarcinoma Of Appendix",
  "Mucinous Adenocarcinoma Of Appendix"
)

## Lung器官的Histologic_category####
histcat_NSCLC <- c(
  "Non-Small Cell Lung Cancer",
  "Non-Small Cell Lung Cancer Recurrent",
  "Non-Small Cell Lung Cancer Metastatic",
  "Her2 Mutant Non-Small Cell Lung Cancer",
  "Large Cell Lung Cancer",
  "Large Cell Lung Cancer Recurrent",
  "Large Cell Lung Cancer Metastatic",
  "Sarcomatoid Carcinoma Of The Lung",
  "Lung Adenocarcinoma Metastatic",
  "Lung Squamous Cell Carcinoma Metastatic",
  "Bronchioloalveolar Carcinoma"
)

histcat_SCLC <- c(
  "Small Cell Lung Cancer",
  "Small Cell Lung Cancer Extensive",
  "Small Cell Lung Cancer Limited",
  "Small Cell Lung Cancer Recurrent",
  "Small Cell Lung Cancer Metastatic",
  "Small Cell Lung Cancer Extensive Stage"
)

## Hematologic器官的Histologic_category####
histcat_Plasma_Cell <- c(
  "Plasma Cell Myeloma",
  "Plasma Cell Myeloma In Remission",
  "Plasma Cell Myeloma Recurrent",
  "Plasma Cell Myeloma Refractory",
  "Plasmacytoma",
  "Plasma Cell Leukaemia",
  "Plasma Cell Leukaemia In Remission",
  "Myeloma Cast Nephropathy"
)

histcat_MDS_MPN <- c(
  "Myelodysplastic Syndrome",
  "Myeloproliferative Neoplasm",
  "Myelodysplastic Syndrome Unclassifiable",
  "Myelodysplastic Syndrome With Excess Blasts",
  "Myelodysplastic Syndrome With Multilineage Dysplasia",
  "Myelodysplastic Syndrome With Single Lineage Dysplasia",
  "Myelodysplastic Syndrome With Ringed Sideroblasts",
  "Myelodysplastic Syndrome Transformation",
  "Myeloproliferative Disorder",
  "Chronic Myelomonocytic Leukaemia",
  "Chronic Myelomonocytic Leukaemia (In Remission)",
  "Juvenile Chronic Myelomonocytic Leukaemia"
)

## Soft tissue器官的Histologic_category####
histcat_GIST <- c(
  "Gastrointestinal Stromal Cancer"
)

## Brain器官的Histologic_category####
histcat_Medulloblastoma <- c(
  "Medulloblastoma",
  "Medulloblastoma Recurrent"
)

## PNS器官的Histologic_category####
histcat_Neuroblastoma <- c(
  "Neuroblastoma",
  "Neuroblastoma Recurrent",
  "Neonatal Neuroblastoma",
  "Ganglioneuroblastoma"
)

## Head/neck器官的Histologic_category####
histcat_Salivary <- c(
  "Salivary Gland Cancer",
  "Her2 Positive Salivary Gland Cancer",
  "Salivary Gland Cancer Recurrent",
  "Metastatic Salivary Gland Cancer",
  "Mucoepidermoid Carcinoma",
  "Adenoid Cystic Carcinoma",
  "Adenoid Cystic Carcinoma Of Salivary Gland",
  "Adenocarcinoma Of Salivary Gland",
  "Acinic Cell Carcinoma Of Salivary Gland",
  "Mucoepidermoid Carcinoma Of Salivary Gland",
  "Carcinoma Ex-Pleomorphic Adenoma",
  "Myoepithelial Carcinoma Of Salivary Gland"
)

## 映射函数：PT → Histologic_category（严格按照用户映射表）####
## 参数说明：
##   - first_pt: PT名称（首字母大写格式）
##   - assigned_type: TCGA癌种代码（可选）
##   - primary_organ: 第一层器官（可选，用于确定默认分类）
## 返回值：用户映射表中定义的Histologic_category
map_to_histologic_category <- function(first_pt, assigned_type = NA_character_, primary_organ = NA_character_) {
  
  # 优先级1：通过TCGA代码映射
  if (!is.na(assigned_type) && assigned_type != "Others" && assigned_type != "") {
    histcat <- tcga_to_histology[[assigned_type]]
    if (!is.null(histcat)) {
      return(histcat)
    }
  }
  
  # 优先级2：通过PT直接映射到Histologic_category
  if (!is.na(first_pt) && first_pt != "") {
    # 统一为Title Case，避免全大写导致的映射失败
    first_pt <- stringr::str_to_title(first_pt)
    # Bowel
    if (first_pt %in% histcat_Small_bowel_cancer) return("Small bowel cancer")
    if (first_pt %in% histcat_Colorectal_cancer) return("Colorectal cancer")
    if (first_pt %in% histcat_Appendiceal_cancer) return("Appendiceal cancer")
    
    # Lung
    if (first_pt %in% histcat_NSCLC) return("NSCLC")
    if (first_pt %in% histcat_SCLC) return("SCLC")
    
    # Hematologic
    if (first_pt %in% histology_Lymphoma) return("Lymphoma")
    if (first_pt %in% histology_Leukemia) return("Leukemia")
    if (first_pt %in% histcat_Plasma_Cell) return("Plasma_Cell")
    if (first_pt %in% histcat_MDS_MPN) return("MDS_MPN")
    
    # Soft tissue
    if (first_pt %in% histology_Sarcoma) return("Soft tissue sarcoma")
    if (first_pt %in% histcat_GIST) return("GIST")
    
    # Brain
    if (first_pt %in% histology_Glioma) return("Glioma")
    if (first_pt %in% histcat_Medulloblastoma) return("Medulloblastoma")
    
    # PNS
    if (first_pt %in% histcat_Neuroblastoma) return("Neuroblastoma")
    
    # Head/neck
    if (first_pt %in% histcat_Salivary) return("Salivary")
    
    # Skin
    if (first_pt %in% histology_Melanoma) return("Melanoma")
  }
  
  # 优先级3：根据器官返回默认的Histologic_category
  if (!is.na(primary_organ) && primary_organ != "") {
    organ_defaults <- list(
      "Head/neck" = "Others",
      "Esophagus/stomach" = "Others",
      "Pancreas" = "Pancreatic cancer",
      "Biliary tract" = "Others",
      "Liver" = "Others",
      "Bowel" = "Others",
      "Breast" = "Breast Cancer",
      "Ovary" = "Ovary cancer",
      "Uterus" = "Others",
      "Cervix" = "Cervical cancer",
      "Prostate" = "Prostate cancer",
      "Kidney" = "Others",
      "Bladder" = "Others",
      "Lung" = "Unknown/Others",
      "Pleura" = "Mesothelioma",
      "Thyroid" = "Thyroid Cancer",
      "Thymus" = "Thymoma",
      "Brain" = "Other",
      "PNS" = "Unknown/Others",
      "Soft tissue" = "Other",
      "Bone" = "Bone cancer",
      "Skin" = "Other",
      "Eye" = "Uveal Melanoma",
      "Hematologic" = "Other",
      "Rare" = "Unknown/Others",
      "CUP" = "Unknown/Others"
    )
    default_cat <- organ_defaults[[primary_organ]]
    if (!is.null(default_cat)) {
      return(default_cat)
    }
  }
  
  # 未匹配到的返回NA
  return(NA_character_)
}

## 仅通过 assigned_type 映射组织学分类（快速版本）####
map_tcga_to_histology <- function(assigned_type) {
  if (is.na(assigned_type) || assigned_type == "Others" || assigned_type == "") {
    return(NA_character_)
  }
  histology <- tcga_to_histology[[assigned_type]]
  if (is.null(histology)) {
    return(NA_character_)
  }
  return(histology)
}

## 仅通过 first_pt 映射组织学分类####
map_pt_to_histology <- function(first_pt) {
  if (is.na(first_pt) || first_pt == "") {
    return(NA_character_)
  }
  # 明确的组织学类型（优先匹配）
  if (first_pt %in% histology_Adenocarcinoma) return("Adenocarcinoma")
  if (first_pt %in% histology_Squamous_Cell_Carcinoma) return("Squamous Cell Carcinoma")
  if (first_pt %in% histology_Urothelial_Carcinoma) return("Urothelial Carcinoma")
  if (first_pt %in% histology_Hepatocellular_Carcinoma) return("Hepatocellular Carcinoma")
  if (first_pt %in% histology_Renal_Cell_Carcinoma) return("Renal Cell Carcinoma")
  if (first_pt %in% histology_Thyroid_Carcinoma) return("Thyroid Carcinoma")
  if (first_pt %in% histology_Basal_Cell_Carcinoma) return("Basal Cell Carcinoma")
  if (first_pt %in% histology_Small_Cell_Carcinoma) return("Small Cell Carcinoma")
  if (first_pt %in% histology_Melanoma) return("Melanoma")
  if (first_pt %in% histology_Sarcoma) return("Sarcoma")
  if (first_pt %in% histology_Glioma) return("Glioma")
  if (first_pt %in% histology_Neuroblastoma) return("Neuroblastoma")
  if (first_pt %in% histology_Lymphoma) return("Lymphoma")
  if (first_pt %in% histology_Leukemia) return("Leukemia")
  if (first_pt %in% histology_Plasma_Cell_Neoplasm) return("Plasma Cell Neoplasm")
  if (first_pt %in% histology_MDS_MPN) return("MDS/MPN")
  if (first_pt %in% histology_Germ_Cell_Tumor) return("Germ Cell Tumor")
  if (first_pt %in% histology_Neuroendocrine_Tumor) return("Neuroendocrine Tumor")
  if (first_pt %in% histology_Embryonal_Tumor) return("Embryonal Tumor")
  if (first_pt %in% histology_Mesothelioma) return("Mesothelioma")
  if (first_pt %in% histology_Thymoma) return("Thymoma")
  if (first_pt %in% histology_Other_Carcinoma) return("Other Carcinoma")
  # Unknown放最后（无法从PT确定组织学类型）
  if (first_pt %in% histology_Unknown) return("Unknown")
  return(NA_character_)
}

## 注意：已删除map_organ_to_histology函数####
## 原设计用于通过primary_organ获取默认组织学分类
## 但所有PT+TCGA代码已完整覆盖第二层映射，不需要此兜底函数

## 获取所有组织学分类列表####
get_histology_list <- function() {
  c(
    # 癌类型
    "Adenocarcinoma",
    "Squamous Cell Carcinoma",
    "Urothelial Carcinoma",
    "Hepatocellular Carcinoma",
    "Renal Cell Carcinoma",
    "Thyroid Carcinoma",
    "Basal Cell Carcinoma",
    "Small Cell Carcinoma",
    "Other Carcinoma",
    # 非癌类型
    "Melanoma",
    "Sarcoma",
    "Glioma",
    "Neuroblastoma",
    "Lymphoma",
    "Leukemia",
    "Plasma Cell Neoplasm",
    "MDS/MPN",
    "Germ Cell Tumor",
    "Neuroendocrine Tumor",
    "Embryonal Tumor",
    "Mesothelioma",
    "Thymoma",
    # 血液系统（泛称）
    "Hematologic Malignancy",
    # 混合/未确定/其他
    "Mixed",
    "Undetermined",
    "Other",
    "Unknown"
  )
}

## 获取组织学分类的大类归属（用于更高层次的分组）####
## 将细分的组织学类型归入更大的类别
get_histology_superclass <- function(histology) {
  carcinoma_types <- c(
    "Adenocarcinoma",
    "Squamous Cell Carcinoma",
    "Urothelial Carcinoma",
    "Hepatocellular Carcinoma",
    "Renal Cell Carcinoma",
    "Thyroid Carcinoma",
    "Basal Cell Carcinoma",
    "Small Cell Carcinoma",
    "Other Carcinoma"
  )
  
  hematologic_types <- c(
    "Lymphoma",
    "Leukemia",
    "Plasma Cell Neoplasm",
    "MDS/MPN",
    "Hematologic Malignancy"
  )
  
  neural_types <- c(
    "Glioma",
    "Neuroblastoma",
    "Neuroendocrine Tumor"
  )
  
  if (is.na(histology) || histology == "") {
    return(NA_character_)
  }
  
  if (histology %in% carcinoma_types) return("Carcinoma")
  if (histology %in% hematologic_types) return("Hematologic")
  if (histology %in% neural_types) return("Neural")
  if (histology == "Melanoma") return("Melanoma")
  if (histology == "Sarcoma") return("Sarcoma")
  if (histology == "Germ Cell Tumor") return("Germ Cell")
  if (histology == "Embryonal Tumor") return("Embryonal")
  if (histology == "Mesothelioma") return("Mesothelioma")
  if (histology == "Thymoma") return("Thymoma")
  if (histology %in% c("Mixed", "Undetermined", "Other", "Unknown")) return("Other/Unknown")
  
  return(NA_character_)
}

## 获取组织学分类映射汇总表####
get_histology_mapping_summary <- function() {
  histology_vars <- c(
    "histology_Adenocarcinoma",
    "histology_Squamous_Cell_Carcinoma",
    "histology_Urothelial_Carcinoma",
    "histology_Hepatocellular_Carcinoma",
    "histology_Renal_Cell_Carcinoma",
    "histology_Thyroid_Carcinoma",
    "histology_Basal_Cell_Carcinoma",
    "histology_Small_Cell_Carcinoma",
    "histology_Other_Carcinoma",
    "histology_Melanoma",
    "histology_Sarcoma",
    "histology_Glioma",
    "histology_Neuroblastoma",
    "histology_Lymphoma",
    "histology_Leukemia",
    "histology_Plasma_Cell_Neoplasm",
    "histology_MDS_MPN",
    "histology_Germ_Cell_Tumor",
    "histology_Neuroendocrine_Tumor",
    "histology_Embryonal_Tumor",
    "histology_Mesothelioma",
    "histology_Thymoma",
    "histology_Unknown"
  )
  
  histology_names <- c(
    "Adenocarcinoma",
    "Squamous Cell Carcinoma",
    "Urothelial Carcinoma",
    "Hepatocellular Carcinoma",
    "Renal Cell Carcinoma",
    "Thyroid Carcinoma",
    "Basal Cell Carcinoma",
    "Small Cell Carcinoma",
    "Other Carcinoma",
    "Melanoma",
    "Sarcoma",
    "Glioma",
    "Neuroblastoma",
    "Lymphoma",
    "Leukemia",
    "Plasma Cell Neoplasm",
    "MDS/MPN",
    "Germ Cell Tumor",
    "Neuroendocrine Tumor",
    "Embryonal Tumor",
    "Mesothelioma",
    "Thymoma",
    "Unknown"
  )
  
  pt_counts <- sapply(histology_vars, function(v) {
    if (exists(v)) length(get(v)) else 0
  })
  
  data.frame(
    histology = histology_names,
    n_pts = pt_counts,
    stringsAsFactors = FALSE
  )
}

# 三层映射对应表数据####
# 生成完整的三层映射表（Primary organs/tissues, Histologic category, TCGA entity）
# 2026-01-19 根据用户校准版本更新
three_layer_mapping <- data.frame(
  Primary_organs_tissues = c(
    # Head/neck（3行）
    "Head/neck", "Head/neck", "Head/neck",
    # Esophagus/stomach（3行）
    "Esophagus/stomach", "Esophagus/stomach", "Esophagus/stomach",
    # Pancreas（1行）
    "Pancreas",
    # Biliary tract（2行）
    "Biliary tract", "Biliary tract",
    # Liver（2行）
    "Liver", "Liver",
    # Bowel（5行）
    "Bowel", "Bowel", "Bowel", "Bowel", "Bowel",
    # Breast（1行）
    "Breast",
    # Ovary（1行）
    "Ovary",
    # Uterus（3行）
    "Uterus", "Uterus", "Uterus",
    # Cervix（1行）
    "Cervix",
    # Prostate（1行）
    "Prostate",
    # Kidney（4行）
    "Kidney", "Kidney", "Kidney", "Kidney",
    # Bladder（2行）
    "Bladder", "Bladder",
    # Lung（4行）
    "Lung", "Lung", "Lung", "Lung",
    # Pleura（1行）
    "Pleura",
    # Thyroid（1行）
    "Thyroid",
    # Thymus（1行）
    "Thymus",
    # Brain（4行）
    "Brain", "Brain", "Brain", "Brain",
    # PNS（3行）
    "PNS", "PNS", "PNS",
    # Soft tissue（3行）
    "Soft tissue", "Soft tissue", "Soft tissue",
    # Bone（1行）
    "Bone",
    # Skin（2行）
    "Skin", "Skin",
    # Eye（1行）
    "Eye",
    # Hematologic（5行）
    "Hematologic", "Hematologic", "Hematologic", "Hematologic", "Hematologic",
    # Rare（3行：ACC, Testicular Cancer, Unknown）
    "Rare", "Rare", "Rare",
    # CUP（1行）
    "CUP"
  ),
  Histologic_category = c(
    # Head/neck
    "Head/neck SCC", "Salivary", "Others",
    # Esophagus/stomach
    "Esophageal Cancer", "Stomach Cancer", "Others",
    # Pancreas
    "Pancreatic cancer",
    # Biliary tract
    "Cholangiocarcinoma", "Others",
    # Liver
    "HCC", "Others",
    # Bowel
    "Small bowel cancer", "Colorectal cancer", "Colorectal cancer", "Appendiceal cancer", "Others",
    # Breast
    "Breast Cancer",
    # Ovary
    "Ovary cancer",
    # Uterus
    "Endometrial", "Sarcoma", "Others",
    # Cervix
    "Cervical cancer",
    # Prostate
    "Prostate cancer",
    # Kidney
    "Renal cell carcinoma", "Renal cell carcinoma", "Renal cell carcinoma", "Others",
    # Bladder
    "Urothelial carcinoma", "Others",
    # Lung
    "NSCLC", "NSCLC", "SCLC", "Unknown/Others",
    # Pleura
    "Mesothelioma",
    # Thyroid
    "Thyroid Cancer",
    # Thymus
    "Thymoma",
    # Brain
    "Glioma", "Glioma", "Medulloblastoma", "Other",
    # PNS
    "PCPG", "Neuroblastoma", "Unknown/Others",
    # Soft tissue
    "Soft tissue sarcoma", "GIST", "Other",
    # Bone
    "Bone cancer",
    # Skin
    "Melanoma", "Other",
    # Eye
    "Uveal Melanoma",
    # Hematologic
    "Lymphoma", "Leukemia", "Plasma_Cell", "MDS_MPN", "Other",
    # Rare
    "ACC", "Testicular Cancer", "Unknown/Others",
    # CUP
    "Unknown/Others"
  ),
  TCGA_entity = c(
    # Head/neck
    "HNSC", NA, NA,
    # Esophagus/stomach
    "ESCA", "STAD", NA,
    # Pancreas
    "PAAD",
    # Biliary tract
    "CHOL", NA,
    # Liver
    "LIHC", NA,
    # Bowel
    NA, "COAD", "READ", NA, NA,
    # Breast
    "BRCA",
    # Ovary
    "OV",
    # Uterus
    "UCEC", "UCS", NA,
    # Cervix
    "CESC",
    # Prostate
    "PRAD",
    # Kidney
    "KIRC", "KIRP", "KICH", NA,
    # Bladder
    "BLCA", NA,
    # Lung
    "LUAD", "LUSC", NA, NA,
    # Pleura
    "MESO",
    # Thyroid
    "THCA",
    # Thymus
    "THYM",
    # Brain
    "GBM", "LGG", NA, NA,
    # PNS
    "PCPG", NA, NA,
    # Soft tissue
    "SARC", NA, NA,
    # Bone
    NA,
    # Skin
    "SKCM", NA,
    # Eye
    "UVM",
    # Hematologic
    "DLBC", "LAML", NA, NA, NA,
    # Rare
    "ACC", "TGCT", NA,
    # CUP
    NA
  ),
  stringsAsFactors = FALSE
)

## 输出三层映射表为CSV####
export_three_layer_mapping <- function(output_path = NULL) {
  if (is.null(output_path)) {
    # 默认输出到Results目录
    output_path <- file.path(dirname(getwd()), "Results", "three_layer_mapping.csv")
  }
  write.csv(three_layer_mapping, output_path, row.names = FALSE, fileEncoding = "UTF-8")
  message(paste("三层映射表已输出:", output_path))
  return(invisible(three_layer_mapping))
}

# 从JADER格式PT中提取英文名称####
# 格式如："日文名(English Name)" 或 "日文名（English Name）" -> "English Name"
# 如果没有括号，保持原样
extract_english_from_jader <- function(pt) {
  if (is.na(pt) || pt == "") {
    return(NA_character_)
  }
  # 尝试匹配半角括号 ()
  match <- regmatches(pt, regexec("\\(([^)]+)\\)$", pt))[[1]]
  if (length(match) > 1 && nchar(match[2]) > 0) {
    # 检查括号内是否是英文（至少包含英文字母）
    if (grepl("[A-Za-z]", match[2])) {
      return(toupper(trimws(match[2])))
    }
  }
  # 尝试匹配全角括号 （）
  match2 <- regmatches(pt, regexec("（([^）]+)）$", pt))[[1]]
  if (length(match2) > 1 && nchar(match2[2]) > 0) {
    # 检查括号内是否是英文（至少包含英文字母）
    if (grepl("[A-Za-z]", match2[2])) {
      return(toupper(trimws(match2[2])))
    }
  }
  # 尝试从字符串中提取任意位置的英文内容（用于复杂格式）
  # 如果PT中有英文+日文混合，优先提取括号内的英文
  # 最后尝试：提取字符串末尾括号内的内容（不限于行末）
  all_matches <- regmatches(pt, gregexpr("[\\(（]([^\\)）]+)[\\)）]", pt, perl = TRUE))[[1]]
  if (length(all_matches) > 0) {
    for (m in rev(all_matches)) {  # 从后往前检查
      inner <- gsub("^[\\(（]|[\\)）]$", "", m)
      if (grepl("[A-Za-z]", inner)) {
        return(toupper(trimws(inner)))
      }
    }
  }
  return(toupper(trimws(pt)))
}

# 智能选择原发癌种PT函数####
# 说明：从多个PT中智能选择代表原发癌种的PT
# 问题：当适应证包含多个PT（如 "METASTASES TO BONE & BREAST CANCER METASTATIC"）时，
#       简单取第一个会选到转移部位PT而非原发癌种PT
# 解决：排除 "Metastases To XXX" 格式的转移部位PT，优先选择原发癌种PT

## 核心函数：从PT列表字符串中选择原发癌种PT####
## 参数说明：
##   - pt_list_str: 用 " & " 连接的PT列表字符串
## 返回值：选定的PT（首字母大写格式，已去除STAGE/GRADE）
## 匹配规则：
##   1. 按 & 分割PT列表
##   2. 排除 "Metastases To XXX" 格式的转移部位PT
##   3. 对剩余PT按优先级排序：能映射到非CUP器官的优先
##   4. 返回优先级最高的PT（标准化后）

select_primary_pt <- function(pt_list_str) {
  if (is.na(pt_list_str) || pt_list_str == "") {
    return(NA_character_)
  }
  
  # 按 & 分割
  pts <- trimws(unlist(strsplit(pt_list_str, "&")))
  
  if (length(pts) == 0) {
    return(NA_character_)
  }
  
  # 对所有PT提取英文名称（转为大写用于后续处理）
  pts <- sapply(pts, extract_english_from_jader, USE.NAMES = FALSE)
  
  # 转移部位PT的模式（Metastases To XXX）
  metastases_to_pattern <- "(?i)^\\s*METASTAS[EI]S\\s+TO\\s+"
  
  # 泛化恶性/转移描述模式（优先级低于具体诊断）
  # 如"Lung Neoplasm Malignant"应让位于"Non-Small Cell Lung Cancer"
  # 如"Metastatic Neoplasm"应让位于"Breast Cancer"
  # 注意：只匹配**完整的**泛化描述，不匹配"器官 CANCER METASTATIC"形式
  # 例如"CANCER METASTATIC"是泛化描述，但"BREAST CANCER METASTATIC"是具体诊断（乳腺癌转移性）
  # 包括：
  #   - NEOPLASM MALIGNANT / MALIGNANT NEOPLASM（完整匹配）
  #   - CANCER METASTATIC / METASTATIC CANCER（完整匹配）
  #   - NEOPLASM METASTATIC / METASTATIC NEOPLASM（完整匹配）
  #   - METASTATIC DISEASE / METASTATIC TUMOUR（完整匹配）
  #   - METASTASIS (单独使用时)
  #   - 器官 + NEOPLASM MALIGNANT（如 LUNG NEOPLASM MALIGNANT，结尾匹配）
  generalized_malignant_pt_pattern <- "(?i)(^\\s*(METASTATIC\\s+(NEOPLASM|CANCER|DISEASE|TUMOU?R)|METASTASIS|NEOPLASM\\s+MALIGNANT|MALIGNANT\\s+NEOPLASM|CANCER\\s+METASTATIC|NEOPLASM\\s+METASTATIC)\\s*$|\\bNEOPLASM\\s+MALIGNANT$)"
  
  # 疾病状态修饰词模式（去除后变为基本癌种名称）
  # 包括：RECURRENT（复发）、REFRACTORY（难治）、IN REMISSION（缓解）
  # LOCALISED/LOCALIZED（局限）、REGIONAL（区域）、EXTENSIVE（广泛）、LIMITED（局限）
  # METASTATIC（转移性）- 去除后保留原发癌种名称
  status_modifiers_pattern <- "(?i)\\s+(RECURRENT|REFRACTORY|IN\\s+REMISSION|LOCALISED|LOCALIZED|REGIONAL|EXTENSIVE\\s+STAGE|EXTENSIVE|LIMITED|\\(IN\\s+REMISSION\\)|METASTATIC)$"
  
  # 标准化函数（去除STAGE/GRADE和状态修饰词，同义词替换，首字母大写）
  normalize_pt <- function(pt) {
    pt_temp <- gsub("(HIGH|LOW)([ -])GRADE", "\\1\\2XPROTECTEDX", pt, ignore.case = TRUE)
    pt_temp <- trimws(sub("\\s*(STAGE|GRADE).*", "", pt_temp, ignore.case = TRUE))
    pt_temp <- gsub("XPROTECTEDX", "GRADE", pt_temp)
    # 去除疾病状态修饰词
    pt_temp <- trimws(sub(status_modifiers_pattern, "", pt_temp, ignore.case = TRUE, perl = TRUE))
    # 同义词替换：将部位特异的癌种名称统一为泛称
    # COLON CANCER / RECTAL CANCER → COLORECTAL CANCER
    pt_temp <- gsub("(?i)^\\s*(COLON|RECTAL)\\s+CANCER\\s*$", "COLORECTAL CANCER", pt_temp, perl = TRUE)
    return(stringr::str_to_title(pt_temp))
  }
  
  # 分离转移部位PT和非转移部位PT
  is_metastases_to <- grepl(metastases_to_pattern, pts, perl = TRUE)
  non_metastases_pts <- pts[!is_metastases_to]
  metastases_pts <- pts[is_metastases_to]
  
  # 如果有非转移部位PT
  if (length(non_metastases_pts) > 0) {
    # 对非转移部位PT进行标准化
    non_metastases_pts_normalized <- sapply(non_metastases_pts, normalize_pt, USE.NAMES = FALSE)
    
    # 识别泛化描述PT
    is_generalized <- grepl(generalized_malignant_pt_pattern, non_metastases_pts, perl = TRUE)
    
    # 去重（如"Plasma Cell Myeloma"和"Plasma Cell Myeloma Recurrent"标准化后相同）
    # 但保留顺序和是否泛化的标记
    unique_mask <- !duplicated(non_metastases_pts_normalized)
    non_metastases_pts_normalized <- non_metastases_pts_normalized[unique_mask]
    is_generalized <- is_generalized[unique_mask]
    
    # 分离具体诊断和泛化描述
    specific_pts <- non_metastases_pts_normalized[!is_generalized]
    generalized_pts <- non_metastases_pts_normalized[is_generalized]
    
    # 优先选择具体诊断中能映射到非CUP器官的PT
    if (length(specific_pts) > 0) {
      for (i in seq_along(specific_pts)) {
        organ <- map_to_primary_organ(specific_pts[i], NA_character_)
        if (!is.na(organ) && organ != "CUP") {
          return(specific_pts[i])
        }
      }
      # 如果具体诊断都是CUP，返回第一个具体诊断
      return(specific_pts[1])
    }
    
    # 如果没有具体诊断，尝试泛化描述
    if (length(generalized_pts) > 0) {
      for (i in seq_along(generalized_pts)) {
        organ <- map_to_primary_organ(generalized_pts[i], NA_character_)
        if (!is.na(organ) && organ != "CUP") {
          return(generalized_pts[i])
        }
      }
      return(generalized_pts[1])
    }
    
    # 如果没有任何非转移部位PT（理论上不会到这里）
    return(non_metastases_pts_normalized[1])
  }
  
  # 如果全是转移部位PT，返回第一个（标准化后）
  return(normalize_pt(pts[1]))
}

## 向量化版本（用于并行处理）####
select_primary_pt_vec <- Vectorize(select_primary_pt)

# irAE SOC（System Organ Class）缩写映射表####
# 说明：为 irAE 列表中的 11 个 SOC 定义更明确、可读性更高的英文缩写
# 用途：用于绘图标签、数据表列名等场景
# 注意：SOC_to_abbrev 用于 SOC全称→缩写，abbrev_to_SOC 用于 缩写→SOC全称

## SOC 全称到缩写的映射####
SOC_to_abbrev <- c(
  "Blood and lymphatic system disorders" = "HEME",    # 血液和淋巴系统疾病
  "Cardiac disorders" = "CARD",                       # 心脏疾病
  "Endocrine disorders" = "ENDO",                     # 内分泌疾病
  "Eye disorders" = "EYE",                            # 眼部疾病
  "Gastrointestinal disorders" = "GI",                # 胃肠道疾病
  "Hepatobiliary disorders" = "HEPA",                 # 肝胆疾病
  "Musculoskeletal and connective tissue disorders" = "MSK",    # 肌肉骨骼和结缔组织疾病
  "Nervous system disorders" = "NEURO",               # 神经系统疾病
  "Renal and urinary disorders" = "RENU",             # 肾脏和泌尿系统疾病
  "Respiratory, thoracic and mediastinal disorders" = "RESP",   # 呼吸、胸部和纵隔疾病
  "Skin and subcutaneous tissue disorders" = "SKIN"   # 皮肤和皮下组织疾病
)

## 缩写到 SOC 全称的反向映射####
abbrev_to_SOC <- setNames(names(SOC_to_abbrev), SOC_to_abbrev)

## 获取 SOC 缩写函数####
## 参数：soc_full - SOC 全称
## 返回：对应的英文缩写，如果未找到则返回 NA
get_soc_abbrev <- function(soc_full) {
  if (is.na(soc_full) || soc_full == "") {
    return(NA_character_)
  }
  abbrev <- SOC_to_abbrev[soc_full]
  if (is.na(abbrev)) {
    return(NA_character_)
  }
  return(as.character(abbrev))
}

## 获取 SOC 全称函数####
## 参数：abbrev - SOC英文缩写
## 返回：对应的 SOC 全称，如果未找到则返回 NA
get_soc_full <- function(abbrev) {
  if (is.na(abbrev) || abbrev == "") {
    return(NA_character_)
  }
  soc_full <- abbrev_to_SOC[abbrev]
  if (is.na(soc_full)) {
    return(NA_character_)
  }
  return(as.character(soc_full))
}

## 获取所有 SOC 缩写列表####
get_soc_abbrev_list <- function() {
  return(SOC_to_abbrev)
}

## 转换数据框列名中的 SOC 为缩写####
## 参数：df - 数据框
## 返回：列名中 SOC_ 前缀后的 SOC 名称被替换为缩写的数据框
rename_soc_cols_to_abbrev <- function(df) {
  old_names <- names(df)
  new_names <- sapply(old_names, function(col_name) {
    # 检查是否是 SOC_ 开头的列
    if (grepl("^SOC_", col_name)) {
      # 提取 SOC 名称（去除 SOC_ 前缀，将下划线替换回空格）
      soc_name <- sub("^SOC_", "", col_name)
      soc_name <- gsub("_", " ", soc_name)
      # 查找缩写
      abbrev <- SOC_to_abbrev[soc_name]
      if (!is.na(abbrev)) {
        return(paste0("SOC_", abbrev))
      }
    }
    return(col_name)
  }, USE.NAMES = FALSE)
  names(df) <- new_names
  return(df)
}

message("映射规则加载完成（包含第一层器官/组织映射、第二层组织学分类映射、三层映射表、智能PT选择函数和irAE SOC缩写映射表）")
