"""
aureonoise — Comprehensive Therapeutic Frequency Map
=====================================================

Research compilation: esoteric, traditional, and non-Western approaches
to therapeutic sound. Every frequency claim with source tradition,
measurement data where available, and evidence quality rating.

Evidence quality key:
    A = peer-reviewed RCT or systematic review
    B = peer-reviewed observational / pilot study
    C = academic / ethnomusicological documentation
    D = traditional system (documented but not empirically tested)
    E = modern practitioner claim (no peer review)

Sources are cited inline. This file is a DATA MODULE — importable as
a Python dict for programmatic access in aureonoise presets, GUI,
and future frequency-selection algorithms.

Compiled: 2026-03-06
"""

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. SOLFEGGIO FREQUENCIES
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Origin: Joseph Puleo (1970s) via Pythagorean numerology on Biblical
# texts. Expanded by Leonard Horowitz in "Healing Codes for the
# Biological Apocalypse" (1999). The 174/285/963 Hz additions are
# post-Horowitz.
#
# Historical claim: Derived from Gregorian chant (Hymn to St. John).
# Reality: Musicologists find no evidence linking specific Hz values
# to medieval solfege. Hertz as a unit did not exist until 1930.
# The original Ut-Re-Mi-Fa-Sol-La were relative pitch names, not
# absolute frequencies.
#
# Numerological pattern: all Solfeggio frequencies reduce to 3, 6, or 9
# via digital root (e.g., 528 -> 5+2+8=15 -> 1+5=6).

SOLFEGGIO = {
    174: {
        "name": "Foundation",
        "solfege": None,  # modern addition, not in original 6
        "claim": "Natural anesthetic; physical/energetic pain relief; security",
        "chakra": None,
        "tradition": "New Age (post-Horowitz expansion)",
        "evidence": "E",
        "digital_root": 3,
    },
    285: {
        "name": "Cellular Repair",
        "solfege": None,
        "claim": "Tissue regeneration, immune system enhancement, cellular memory",
        "chakra": None,
        "tradition": "New Age (post-Horowitz expansion)",
        "evidence": "E",
        "digital_root": 6,
    },
    396: {
        "name": "Liberating Guilt and Fear",
        "solfege": "UT",
        "claim": "Turning grief into joy; releasing guilt and fear",
        "chakra": "Root (Muladhara)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "digital_root": 9,
    },
    417: {
        "name": "Undoing Situations and Facilitating Change",
        "solfege": "RE",
        "claim": "Clearing negative memories; enabling perspective shifts",
        "chakra": "Sacral (Svadhisthana)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "digital_root": 3,
    },
    528: {
        "name": "Transformation and Miracles (Love Frequency)",
        "solfege": "MI",
        "claim": "DNA repair, increased energy, clarity, creativity",
        "chakra": "Solar Plexus (Manipura)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "note": (
            "One in-vitro study showed 528 Hz reduced ethanol-induced "
            "cell death in astrocytes. This does NOT demonstrate DNA repair. "
            "No peer-reviewed evidence for direct DNA repair exists."
        ),
        "digital_root": 6,
    },
    639: {
        "name": "Reconnecting and Balancing Relationships",
        "solfege": "FA",
        "claim": "Social harmony, communication, understanding",
        "chakra": "Heart (Anahata)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "digital_root": 9,
    },
    741: {
        "name": "Solving Problems and Self-Expression",
        "solfege": "SOL",
        "claim": "Emotional expression, toxin removal, intuition",
        "chakra": "Throat (Vishuddha)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "digital_root": 3,
    },
    852: {
        "name": "Awakening Intuition",
        "solfege": "LA",
        "claim": "Returning to spiritual order; higher consciousness",
        "chakra": "Third Eye (Ajna)",
        "tradition": "Puleo/Horowitz (1999)",
        "evidence": "E",
        "digital_root": 6,
    },
    963: {
        "name": "Higher Consciousness / God Frequency",
        "solfege": None,
        "claim": "Crown chakra activation, pineal gland, 'divine connection'",
        "chakra": "Crown (Sahasrara)",
        "tradition": "New Age (post-Horowitz expansion)",
        "evidence": "E",
        "digital_root": 9,
    },
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. PYTHAGOREAN TUNING & HEALING INTERVALS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Pythagoras (c. 570-495 BCE) discovered that consonant intervals
# correspond to simple integer ratios via monochord experiments.
# The "healing" application is a modern overlay.

PYTHAGOREAN_INTERVALS = {
    "unison": {
        "ratio": (1, 1),
        "cents": 0,
        "claim": "Unity, grounding",
        "tradition": "Greek (Pythagorean)",
        "evidence": "C",
    },
    "octave": {
        "ratio": (2, 1),
        "cents": 1200,
        "claim": "Completeness, return to source",
        "tradition": "Greek (Pythagorean)",
        "evidence": "C",
    },
    "perfect_fifth": {
        "ratio": (3, 2),
        "cents": 701.96,
        "claim": "Most consonant interval after octave; 'soul alignment'; "
                 "used in C256/G384 tuning forks for nervous system balance",
        "tuning_fork_pair": (256, 384),
        "tradition": "Greek (Pythagorean) / modern sound healing",
        "evidence": "C/E",
    },
    "perfect_fourth": {
        "ratio": (4, 3),
        "cents": 498.04,
        "claim": "Stability, foundation",
        "tradition": "Greek (Pythagorean)",
        "evidence": "C",
    },
    "major_third": {
        "ratio": (5, 4),
        "cents": 386.31,
        "claim": "Joy, warmth (just intonation, not strict Pythagorean)",
        "tradition": "Renaissance extension of Pythagorean",
        "evidence": "C",
    },
    "minor_third": {
        "ratio": (6, 5),
        "cents": 315.64,
        "claim": "Melancholy, introspection",
        "tradition": "Renaissance extension",
        "evidence": "C",
    },
}

# Pythagorean tuning fork base frequencies used in sound healing
PYTHAGOREAN_TUNING_FORKS = {
    "C": 256.0,   # "Scientific pitch" / Verdi C
    "D": 288.0,
    "E": 320.0,
    "F": 341.3,
    "G": 384.0,
    "A": 426.7,
    "B": 480.0,
    "C_high": 512.0,
    "note": "Based on C4=256 Hz (scientific/philosophical pitch), NOT A440 standard",
    "tradition": "Neo-Pythagorean sound healing",
    "evidence": "E",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. INDIAN RAGA THERAPY (RAAG CHIKITSA)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Indian classical music uses relative pitch (Sa can be any frequency).
# The system works through melodic contour (raga), not absolute Hz.
# 72 melakarta ragas (Carnatic) are claimed to correspond to 72 nerves.
#
# Swara ratios (22-shruti system, Sa=1):
#   Sa=1/1, Ri(komal)=256/243, Ri=9/8, Ga(komal)=32/27, Ga=5/4,
#   Ma=4/3, Ma(tivra)=45/32, Pa=3/2, Dha(komal)=128/81, Dha=5/3,
#   Ni(komal)=16/9, Ni=15/8
#
# Time Theory (Samay Chakra): ragas assigned to 8 prahars (3-hr periods)
# based on the vadi (dominant) swara and madhyam type.

SWARA_RATIOS = {
    "Sa":         {"ratio": (1, 1),     "cents": 0,      "western": "C (tonic)"},
    "Ri_komal":   {"ratio": (256, 243), "cents": 90.22,  "western": "Db"},
    "Ri_shuddha": {"ratio": (9, 8),     "cents": 203.91, "western": "D"},
    "Ga_komal":   {"ratio": (32, 27),   "cents": 294.13, "western": "Eb"},
    "Ga_shuddha": {"ratio": (5, 4),     "cents": 386.31, "western": "E"},
    "Ma_shuddha": {"ratio": (4, 3),     "cents": 498.04, "western": "F"},
    "Ma_tivra":   {"ratio": (45, 32),   "cents": 590.22, "western": "F#"},
    "Pa":         {"ratio": (3, 2),     "cents": 701.96, "western": "G"},
    "Dha_komal":  {"ratio": (128, 81),  "cents": 792.18, "western": "Ab"},
    "Dha_shuddha":{"ratio": (5, 3),     "cents": 884.36, "western": "A"},
    "Ni_komal":   {"ratio": (16, 9),    "cents": 996.09, "western": "Bb"},
    "Ni_shuddha": {"ratio": (15, 8),    "cents": 1088.27,"western": "B"},
}

# If Sa=240 Hz (common reference):
SWARA_HZ_SA240 = {
    "Sa": 240.0, "Ri_komal": 252.8, "Ri_shuddha": 270.0,
    "Ga_komal": 284.4, "Ga_shuddha": 300.0, "Ma_shuddha": 320.0,
    "Ma_tivra": 337.5, "Pa": 360.0, "Dha_komal": 379.3,
    "Dha_shuddha": 400.0, "Ni_komal": 426.7, "Ni_shuddha": 450.0,
}

RAGA_THERAPY = {
    # Format: raga_name -> {time, condition, emotion, notes}
    # Source: ICMACY Raga Chikitsa compendium + Manasukh Dhvani
    "Ahir Bhairav":     {"time": "morning", "conditions": ["indigestion", "rheumatic arthritis", "hypertension"], "emotion": "compassion"},
    "Asavari":          {"time": "morning", "conditions": [], "emotion": "confidence"},
    "Bageshri":         {"time": "late night", "conditions": ["insomnia"], "emotion": "rest"},
    "Bhairavi":         {"time": "early morning", "conditions": ["rheumatoid arthritis", "sinusitis"], "emotion": "peace, celebration, detachment"},
    "Bhimpalasi":       {"time": "afternoon", "conditions": ["anxiety", "hypertension"], "emotion": "success"},
    "Bhupali":          {"time": "evening", "conditions": ["tension", "anger", "mental fatigue"], "emotion": "peace"},
    "Brindavani Sarang":{"time": "midday", "conditions": ["depression"], "emotion": "wisdom, energy"},
    "Chandrakauns":     {"time": "night", "conditions": ["anorexia"], "emotion": "normalize weight"},
    "Darbari Kanada":   {"time": "late night", "conditions": ["headache", "asthma"], "emotion": "calmness, mental ease"},
    "Durga":            {"time": "evening", "conditions": [], "emotion": "joy, compassion, self-confidence"},
    "Gunakali":         {"time": "morning", "conditions": ["rheumatic arthritis", "constipation", "headache", "hemorrhoids"], "emotion": "settled mind"},
    "Hansadhwani":      {"time": "evening", "conditions": [], "emotion": "celebration, happiness"},
    "Hindol":           {"time": "morning", "conditions": ["rheumatic arthritis", "spondylitis", "backache", "hypertension"], "emotion": "devotion"},
    "Jaunpuri":         {"time": "afternoon", "conditions": ["intestinal gas", "diarrhea", "constipation"], "emotion": "relief"},
    "Kafi":             {"time": "night", "conditions": ["sleep disorders"], "emotion": "creativity"},
    "Kalyan":           {"time": "evening", "conditions": ["joint maintenance"], "emotion": "compassion"},
    "Kedar":            {"time": "night", "conditions": ["headache", "cold", "cough", "asthma", "sleep disorders"], "emotion": "devotion"},
    "Khamaj":           {"time": "night", "conditions": ["headache", "sleep disorders"], "emotion": "calmness"},
    "Malhar":           {"time": "rainy season", "conditions": ["asthma", "sunstroke"], "emotion": "joy"},
    "Malkauns":         {"time": "late night", "conditions": ["intestinal gas"], "emotion": "tranquility, restful sleep"},
    "Marwa":            {"time": "evening", "conditions": ["indigestion", "hyperacidity"], "emotion": "coherence"},
    "Puriya":           {"time": "evening", "conditions": ["colitis", "anemia", "hypertension"], "emotion": "harmony, peace"},
    "Puriya Dhanashri": {"time": "evening", "conditions": ["anemia"], "emotion": "relaxation"},
    "Rageshri":         {"time": "night", "conditions": [], "emotion": "rejuvenation, longevity"},
    "Shree":            {"time": "evening", "conditions": ["anorexia", "cold", "cough", "asthma"], "emotion": "devotion"},
    "Todi":             {"time": "morning", "conditions": ["blood pressure normalization"], "emotion": "joy"},
    "Yaman":            {"time": "early night", "conditions": ["rheumatic arthritis"], "emotion": "compassion, joy"},
}

RAGA_THERAPY_META = {
    "tradition": "Ayurveda / Nada Yoga / Indian Classical (Hindustani + Carnatic)",
    "evidence": "B/C",
    "note": (
        "Modern studies show Indian classical music reduces cortisol, "
        "activates parasympathetic nervous system, and induces alpha/theta "
        "brainwave states. Recent research suggests efficacy in Alzheimer's, "
        "dementia, and autism. Neural firing increases during raga listening. "
        "72 melakarta ragas -> 72 nerve claim is traditional, not empirically verified."
    ),
    "key_principle": "Melodic contour (raga) drives the effect, not absolute frequency",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. CHINESE FIVE ELEMENTS SOUND HEALING (WU XING / FPMT)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Five Phases Music Therapy (FPMT) maps pentatonic notes to the
# Wu Xing (Five Elements), organs, emotions, seasons, and colors.
#
# The character for 'medicine' (藥 yao) derives from 'music' (樂 yue).
#
# Therapeutic principle: emotion is subdued by the conquering element.
# e.g., Anger (Wood) subdued by Metal (Shang tone).

CHINESE_FIVE_ELEMENTS = {
    "Gong": {
        "pinyin": "Gōng (宫)",
        "element": "Earth (土)",
        "organ": "Spleen",
        "emotion_positive": "Thinking/Reflection",
        "emotion_excess": "Worry/Overthinking",
        "season": "Late Summer",
        "western_key": "C",
        "approximate_hz": 256,  # C4 in scientific pitch
        "color": "Yellow",
        "direction": "Center",
        "taste": "Sweet",
    },
    "Shang": {
        "pinyin": "Shāng (商)",
        "element": "Metal (金)",
        "organ": "Lungs",
        "emotion_positive": "Courage",
        "emotion_excess": "Grief/Sorrow",
        "season": "Autumn",
        "western_key": "D",
        "approximate_hz": 288,
        "color": "White",
        "direction": "West",
        "taste": "Pungent",
    },
    "Jue": {
        "pinyin": "Jué (角)",
        "element": "Wood (木)",
        "organ": "Liver",
        "emotion_positive": "Benevolence",
        "emotion_excess": "Anger",
        "season": "Spring",
        "western_key": "E",
        "approximate_hz": 320,
        "color": "Green/Blue-Green",
        "direction": "East",
        "taste": "Sour",
    },
    "Zhi": {
        "pinyin": "Zhǐ (徵)",
        "element": "Fire (火)",
        "organ": "Heart",
        "emotion_positive": "Joy/Propriety",
        "emotion_excess": "Overexcitement/Mania",
        "season": "Summer",
        "western_key": "G",
        "approximate_hz": 384,
        "color": "Red",
        "direction": "South",
        "taste": "Bitter",
    },
    "Yu": {
        "pinyin": "Yǔ (羽)",
        "element": "Water (水)",
        "organ": "Kidneys",
        "emotion_positive": "Wisdom",
        "emotion_excess": "Fear/Fright",
        "season": "Winter",
        "western_key": "A",
        "approximate_hz": 427,
        "color": "Black/Dark Blue",
        "direction": "North",
        "taste": "Salty",
    },
}

CHINESE_FIVE_ELEMENTS_META = {
    "tradition": "Traditional Chinese Medicine (TCM) / Wu Xing",
    "evidence": "B/C",
    "therapeutic_principle": (
        "Five-Element conquest cycle subdues excess emotion: "
        "Wood conquers Earth, Earth conquers Water, Water conquers Fire, "
        "Fire conquers Metal, Metal conquers Wood. "
        "E.g., excess anger (Wood) -> listen to Shang (Metal)."
    ),
    "note": (
        "Research in FPMT (Five Phases Music Therapy) shows measurable "
        "effects on stress, insomnia, and emotional regulation. "
        "PMC articles document AI-assisted diagnostics using five-tone analysis."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. TIBETAN SINGING BOWLS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Hand-hammered metal bowls (typically 7-metal alloy).
# Key acoustic feature: inherent binaural beats from a single source.
# The two vibrational modes of the bowl wall produce slightly different
# frequencies, creating MONOPHONIC BINAURAL BEATS.
#
# Characteristic interval: flatted fifth (tritone / augmented fourth).

TIBETAN_SINGING_BOWLS = {
    "frequency_range": (100, 900),  # Hz, typical range
    "fundamental_range": (100, 500),  # Hz, fundamentals only
    "characteristic_interval": "tritone (flatted fifth / augmented fourth)",
    "harmonic_series": {
        "typical_partials": [3, 6, 10, 14],  # harmonics present in 6-12" bowls
        "extended_partials": [2, 3, 4, 5, 6, 9, 10, 12, 14],  # some bowls
    },
    "binaural_beat_mechanism": (
        "Each vibrational mode of the bowl wall splits into two close "
        "frequencies due to slight asymmetry. The interference creates "
        "pulsation (beating). This is monophonic binaural beating — "
        "both frequencies from one source, not two separated sources."
    ),
    "brainwave_entrainment": {
        "beta_to_alpha": "Initial state transition",
        "alpha_to_theta": "Deep meditation, some bowls",
    },
    "measurement_challenge": (
        "Frequencies are 'moving targets' — readings spread over 1-5 Hz range. "
        "Software reports average or mode of hundreds of readings/second."
    ),
}

SINGING_BOWL_CHAKRA_MAP = {
    # Standard (approximate) assignments used in sound healing
    "C": {"hz_approx": 256, "chakra": "Root (Muladhara)",      "claim": "Grounding, stability"},
    "D": {"hz_approx": 288, "chakra": "Sacral (Svadhisthana)", "claim": "Creativity, sexuality"},
    "E": {"hz_approx": 320, "chakra": "Solar Plexus (Manipura)","claim": "Willpower, confidence"},
    "F": {"hz_approx": 341, "chakra": "Heart (Anahata)",       "claim": "Love, compassion"},
    "G": {"hz_approx": 384, "chakra": "Throat (Vishuddha)",    "claim": "Communication, expression"},
    "A": {"hz_approx": 427, "chakra": "Third Eye (Ajna)",      "claim": "Intuition, insight"},
    "B": {"hz_approx": 480, "chakra": "Crown (Sahasrara)",     "claim": "Spiritual connection"},
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. DIDGERIDOO THERAPY
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Australian Aboriginal instrument (yidaki). Therapeutic mechanism is
# UPPER AIRWAY MUSCLE TRAINING, not acoustic frequency per se.
# Landmark RCT: Puhan et al., BMJ 2006.

DIDGERIDOO = {
    "fundamental_hz": {
        "range": (60, 90),
        "typical": 70,
        "note_range": "B1 to F2",
        "measured_instruments": {
            # From UNSW Physics (Joe Wolfe lab)
            "arnhem_land_1": 60,    # B1
            "arnhem_land_2": 80,    # E2
            "arnhem_land_3": 64,    # C2
        },
    },
    "formant_range": (500, 3500),  # Hz, shaped by mouth cavity
    "resonance_range": (1000, 2000),  # Hz, important spectral shaping
    "harmonic_content": "All harmonics of fundamental; ~12 dB/octave rolloff",
    "playing_pressure": {
        "drone_kpa": (1, 2),
        "second_mode_kpa": (4, 5),
    },
    "mode_frequency_ratio": {
        "second_to_first": (1.30, 1.43),  # vs 1.50 for perfect cylinder
    },
    "sleep_apnea_rct": {
        "source": "Puhan et al., BMJ 2006; doi:10.1136/bmj.38705.470590.55",
        "n": 25,
        "intervention_group": 14,
        "control_group": 11,
        "duration_months": 4,
        "practice_days_per_week": 5.9,
        "practice_minutes_per_day": 25.3,
        "instrument": "Standardized acrylic didgeridoo, 130 cm, 4 cm diameter",
        "primary_outcome_ess": {
            "didgeridoo_change": -4.4,
            "control_change": -1.4,
            "difference": -3.0,
            "ci_95": (-5.7, -0.3),
            "p": 0.03,
        },
        "ahi_change": {
            "didgeridoo": -10.7,
            "control": -4.5,
            "difference": -6.2,
            "ci_95": (-12.3, -0.1),
            "p": 0.05,
        },
        "partner_sleep_disturbance": {
            "difference": -2.8,
            "ci_95": (-4.7, -0.9),
            "p_lt": 0.01,
        },
        "evidence": "A",
    },
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. SUFI SOUND HEALING (DHIKR)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Dhikr ("remembrance") = repetitive recitation of divine names/phrases.
# No specific Hz values documented in the tradition — the system works
# through phonemic vibration, breath rhythm, and attentional focus.
# Research shows alpha wave induction during dhikr.

SUFI_DHIKR = {
    "primary_sounds": {
        "Hu": {
            "pronunciation": "whoo",
            "body_focus": "Solar plexus",
            "claim": "Source of life/light/love; fills body with vital energy",
            "vowel": "U (close back rounded)",
        },
        "La_ilaha_illa_Allah": {
            "claim": "Complete body/mind purification through rhythmic repetition",
            "breath_pattern": "Exhale-inhale cycle synchronized with phrase",
        },
        "Ya_Hayy": {
            "claim": "Invocation of 'The Living'; vitality, life force",
        },
        "Ya_Shafi": {
            "claim": "Invocation of 'The Healer'; specific healing intent",
        },
    },
    "measured_effects": {
        "brainwave": "Alpha wave induction via auditory nerve stimulation",
        "mechanism": (
            "Dhikrullah through the auditory nerve stimulates the brain "
            "to present alpha waves. Strong relationship between dhikr "
            "vibration and brain wave patterns documented."
        ),
        "source": "Bircu Journal; Scitepress 2019",
    },
    "tradition": "Islamic Sufism (Ottoman Empire musical tradition particularly)",
    "evidence": "B/C",
    "note": (
        "No spectrographic frequency measurements of specific dhikr sounds "
        "found in peer-reviewed literature. The system operates through "
        "phonemic resonance and attentional/devotional focus, not Hz targeting."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 8. SHAMANIC DRUMMING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Rhythmic percussion for trance induction. The key parameter is
# BEAT RATE (repetition frequency), not the acoustic frequency of
# the drum itself.

SHAMANIC_DRUMMING = {
    "beat_rate_hz": {
        "range": (3.0, 7.0),
        "optimal": 4.5,
        "harner_technique": 3.67,  # 220 BPM / 60 = 3.67 Hz
        "unit": "beats per second (Hz)",
    },
    "brainwave_target": {
        "band": "Theta",
        "range_hz": (4.0, 7.0),
        "peak_response": 4.5,
    },
    "eeg_research": {
        "maxfield_1990": {
            "finding": "4.5 bps drumming -> strong theta increase + alpha + beta",
            "n": 12,
            "evidence": "B",
        },
        "neher_1962": {
            "finding": "Auditory driving: brain frequencies synchronize with "
                       "rhythmic drumming at 3-8 Hz",
            "n": 10,
            "evidence": "B",
        },
        "flor_henry_2017": {
            "finding": "Shamanic trance shows altered hemispheric laterality, "
                       "systemic psychobiology changes on high-density EEG",
            "source": "Cogent Psychology 2017",
            "evidence": "B",
        },
        "hove_2016": {
            "finding": "Neural correlates of shamanic state of consciousness "
                       "documented via EEG",
            "source": "PMC 8012721",
            "evidence": "B",
        },
    },
    "drum_acoustic_frequency": {
        "note": "The drum's acoustic pitch (typically 50-200 Hz fundamental) "
                "is secondary. The therapeutic variable is the REPETITION RATE.",
    },
    "tradition": "Cross-cultural shamanic (Siberian, Native American, Nordic, Amazonian)",
    "evidence": "B",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 9. CRYSTAL SINGING BOWLS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Made from crushed quartz crystal (99.99% SiO2). Produce purer,
# more sustained tones than metal bowls. Often tuned to specific
# pitches at the factory.

CRYSTAL_SINGING_BOWLS = {
    "material": "Quartz crystal (SiO2, 99.99% purity)",
    "frequency_range": (100, 900),  # Hz
    "characteristics": [
        "Purer fundamental than metal bowls (fewer inharmonic partials)",
        "Longer sustain",
        "Precise factory tuning possible",
        "Available in 432 Hz and 440 Hz tuning systems",
    ],
    "tuning_systems": {
        "440_hz": {"A4": 440.0, "note": "Standard concert pitch"},
        "432_hz": {"A4": 432.0, "note": "Claimed 'natural' tuning; Verdi's preference"},
    },
    "solfeggio_sets": {
        "description": "Sets tuned to Solfeggio frequencies (396-852 Hz)",
        "frequencies": [396, 417, 528, 639, 741, 852],
    },
    "chakra_mapping": "Same as SINGING_BOWL_CHAKRA_MAP above",
    "tradition": "Modern sound healing (late 20th century invention)",
    "evidence": "E",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 10. PLANETARY FREQUENCIES (HANS COUSTO — COSMIC OCTAVE)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Hans Cousto (Swiss mathematician, 1978): orbital periods of celestial
# bodies octave-transposed into the audible range using f * 2^n.
#
# The math is deterministic and reproducible. The therapeutic claims
# layered on top are not empirically validated.

PLANETARY_FREQUENCIES = {
    "Sun": {
        "hz": 126.22,
        "note": "C (approx)",
        "color": "Green to yellow-green",
        "wavelength_nm": 540,
        "claim": "Centering, magic, transcendental awareness",
        "chakra": None,
    },
    "Earth_day": {
        "hz": 194.18,
        "note": "G",
        "color": "Orange-red",
        "wavelength_nm": 700,
        "chakra": "Root (Muladhara) / Kundalini",
        "claim": "Dynamic, vitalizing; tonifying (medicinal)",
    },
    "Earth_year": {
        "hz": 136.10,
        "note": "C#",
        "color": "Blue-green (turquoise)",
        "wavelength_nm": 500,
        "chakra": "Heart (Anahata)",
        "claim": "Relaxing, soothing, balancing; sedative (medicinal)",
        "note_extra": (
            "This frequency is identical to the fundamental of Indian "
            "classical music (Sa) in the tradition where Sa=136.1 Hz, "
            "and to the Tibetan OM. Corresponds to A=432.10 Hz concert pitch."
        ),
    },
    "Earth_platonic_year": {
        "hz": 172.06,
        "note": "F",
        "color": "Purple/Violet",
        "wavelength_nm": 400,
        "chakra": "Crown (Sahasrara)",
        "claim": "Cheerfulness, clarity of spirit; anti-depressive",
    },
    "Moon_synodic": {
        "hz": 210.42,
        "note": "G#",
        "color": "Orange",
        "wavelength_nm": 650,
        "chakra": "Sacral (Svadhisthana)",
        "claim": "Sexual energy, erotic communication; menstruation regulation",
    },
    "Moon_sidereal": {
        "hz": 227.43,
        "note": "A#",
        "color": "Yellow",
    },
    "Moon_culmination": {
        "hz": 187.61,
        "note": "F#",
        "color": "Red",
    },
    "Mercury": {
        "hz": 141.27,
        "note": "C#/D",
        "color": "Blue-green",
        "wavelength_nm": 480,
        "chakra": "Throat (Vishuddha)",
        "claim": "Speech center, communicative-intellectual principle",
    },
    "Venus": {
        "hz": 221.23,
        "note": "A",
        "color": "Yellow-orange",
        "wavelength_nm": 615,
        "chakra": "Third Eye (Ajna)",
        "claim": "Higher love energy, aspiration for harmony",
    },
    "Mars": {
        "hz": 144.72,
        "note": "D",
        "color": "Blue",
        "wavelength_nm": 470,
        "claim": "Strength of will, focused energy",
    },
    "Jupiter": {
        "hz": 183.58,
        "note": "F#",
        "color": "Red",
        "wavelength_nm": 740,
        "claim": "Creative power, continuous construction",
    },
    "Saturn": {
        "hz": 147.85,
        "note": "D",
        "color": "Blue",
        "wavelength_nm": 460,
        "claim": "Concentration, becoming conscious; 'cosmic controller'",
    },
    "Uranus": {
        "hz": 207.36,
        "note": "G#",
        "color": "Orange",
        "wavelength_nm": 560,
        "claim": "Surprise, renewal, primeval/erotic power",
    },
    "Neptune": {
        "hz": 211.44,
        "note": "G#",
        "color": "Orange",
        "wavelength_nm": 645,
        "claim": "Intuition, unconscious, dream experience",
    },
    "Pluto": {
        "hz": 140.25,
        "note": "C#",
        "color": "Blue-green",
        "wavelength_nm": 485,
        "claim": "Group dynamic, integration into societal structures",
    },
    "Earth_sidereal_day": {
        "hz": 194.71,
        "note": "G",
        "color": "Red-orange",
    },
    # Extended bodies (from product catalogs)
    "Sedna":    {"hz": 128.10},
    "Chiron":   {"hz": 151.27},
    "Nibiru":   {"hz": 161.26},  # speculative
}

PLANETARY_META = {
    "author": "Hans Cousto",
    "year": 1978,
    "formula": "f_audible = f_orbital * 2^n (octave transposition into 20-20000 Hz)",
    "tradition": "Western esoteric / mathematical cosmology",
    "evidence": "C (math is valid; therapeutic claims are E)",
    "book": "The Cosmic Octave: Origin of Harmony (1978)",
    "note": (
        "The calculation method is deterministic: given an orbital period, "
        "the resulting frequency is mathematically fixed. The assignment of "
        "therapeutic/psychological properties to these frequencies is "
        "Cousto's interpretive layer, not an empirical finding."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 11. 528 Hz "DNA REPAIR" — EVIDENCE ASSESSMENT
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DNA_528_CLAIMS = {
    "frequency_hz": 528,
    "popular_name": "Love Frequency / Miracle Tone / DNA Repair Frequency",
    "claim": "Repairs DNA, promotes healing at cellular level",
    "origin": "Leonard Horowitz, 'Healing Codes for the Biological Apocalypse' (1999)",
    "evidence_assessment": {
        "in_vitro_astrocyte": {
            "finding": "528 Hz exposure reduced ethanol-induced cell death in astrocytes",
            "note": "Cell death reduction != DNA repair. Specific to one cell type in vitro.",
            "evidence": "B (limited, preliminary)",
        },
        "dna_repair_direct": {
            "finding": "No peer-reviewed study demonstrates direct DNA repair from 528 Hz",
            "evidence": "No evidence",
        },
        "relaxation_anxiety": {
            "finding": "Some studies show solfeggio tones reduce anxiety and promote relaxation",
            "note": "Non-specific to 528 Hz; could be any calming sound",
            "evidence": "B (weak)",
        },
    },
    "scientific_consensus": (
        "The 'DNA repair' claim is not supported by scientific evidence. "
        "DNA does not have a resonant frequency at 528 Hz. The term 'repair' "
        "implies a specific molecular mechanism that has not been demonstrated."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 12. SCHUMANN RESONANCE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Electromagnetic resonances of the Earth-ionosphere cavity.
# Discovered by Winfried Otto Schumann (1952).
# These are REAL, measurable electromagnetic phenomena.

SCHUMANN_RESONANCE = {
    "harmonics_hz": [7.83, 14.3, 20.8, 27.3, 33.8],
    "fundamental": 7.83,
    "description": "ELF electromagnetic resonances in Earth-ionosphere cavity",
    "discovery": "W.O. Schumann, 1952 (predicted) / 1954 (measured)",
    "brainwave_correspondence": {
        7.83:  "Alpha/Theta border — relaxed wakefulness, light meditation",
        14.3:  "Low Beta — alert, focused consciousness",
        20.8:  "Beta — active thinking",
        27.3:  "High Beta — high arousal",
        33.8:  "Gamma border — cognitive processing",
    },
    "eeg_band_overlap": {
        "delta": (0.5, 4.0),
        "theta": (4.0, 8.0),
        "alpha": (8.0, 13.0),
        "beta": (13.0, 30.0),
        "gamma": (30.0, 100.0),
        "note": "Schumann harmonics span the same range as the first 4 EEG bands",
    },
    "health_claims": {
        "circadian_rhythm": "May help stabilize circadian rhythms",
        "anxiety_reduction": "Exposure to 7.83 Hz fields may reduce anxiety",
        "memory_enhancement": "Some evidence for improved memory",
        "cellular_mechanism": (
            "ELF may modulate cellular calcium influx/efflux via "
            "field-sensitive molecules or radical pairs affecting ion channels"
        ),
    },
    "tradition": "Western physics / biophysics",
    "evidence": "B (well-documented physical phenomenon; health correlations preliminary)",
    "note": (
        "The Schumann resonances are real EM phenomena, measurable with "
        "magnetometers. The correlation with EEG bands is a frequency overlap, "
        "not a proven causal mechanism. Further research needed."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 13. SACRED GEOMETRY & GOLDEN RATIO IN SOUND
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Phi (φ = 1.618033988749895) appears in music as interval ratios,
# Fibonacci-derived rhythmic structures, and spectral relationships.

PHI_SOUND = {
    "phi": 1.618033988749895,
    "phi_as_frequency": {
        "hz": 1.618,
        "claim": "Golden Ratio meditation frequency; autonomic nervous system balance",
        "evidence": "E",
        "note": "1.618 Hz is below audible range; used as binaural/isochronic beat rate",
    },
    "phi_interval": {
        "cents": 833.09,  # 1200 * log2(phi)
        "nearest_interval": "Between minor sixth (814c) and major sixth (884c)",
        "claim": "Most aesthetically balanced interval; bridges consonance/dissonance",
    },
    "solfeggio_phi_relationships": {
        "description": (
            "Multiplying consecutive Solfeggio frequencies by phi yields "
            "approximate relationships: 396 * 1.618 ≈ 640.7 (close to 639), "
            "528 * 1.618 ≈ 854.3 (close to 852)"
        ),
    },
    "fibonacci_frequencies": {
        "sequence": [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597],
        "claim": "Fibonacci Hz values carry phi-proportioned energy",
        "evidence": "E",
    },
    "A432_phi_connection": {
        "claim": "432 Hz = C5 (512 Hz) / phi^(some power) — various numerological derivations",
        "reality": "432 / 512 = 0.84375 ≠ any power of phi. The connection is forced.",
    },
    "tradition": "Western esoteric / sacred geometry",
    "evidence": "E (mathematical relationships exist; therapeutic claims unsubstantiated)",
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 14. BIOACOUSTIC BIOLOGY (Sharry Edwards)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Voice spectral analysis → identify "missing" frequencies →
# deliver complementary low-frequency sound to restore balance.
# Pioneer: Sharry Edwards, Institute of BioAcoustic Biology.

BIOACOUSTIC_BIOLOGY = {
    "method": (
        "Voice spectral analysis identifies frequency deficits. "
        "Complementary low-frequency sounds are delivered back to the body "
        "to 'harmonize' the imbalance."
    ),
    "claims": [
        "Control pain via specific frequencies",
        "Regulate body temperature",
        "Normalize heart rhythm and blood pressure",
        "Regenerate body tissue",
        "Address macular degeneration, MS, headaches, brain trauma",
    ],
    "pioneer": "Sharry Edwards, M.Ed.",
    "institution": "Institute of BioAcoustic Biology, Albany, OH",
    "recognition": "Included in Duke University Encyclopedia of New Medicine",
    "frequency_range": "Low-frequency sound (specific Hz not publicly documented)",
    "tradition": "Modern Western alternative medicine",
    "evidence": "D/E",
    "note": (
        "No peer-reviewed RCTs found. Claims are broad and extraordinary. "
        "The voice analysis methodology and treatment protocols are "
        "proprietary, limiting independent replication."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 15. ORGAN RESONANCE FREQUENCIES (Barbara Hero / Lambdoma)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Barbara Hero calculated resonant frequencies for human organs using
# the Lambdoma matrix (Pythagorean ratio grid based on 256 Hz).
# Method: sound passed through organs + mathematical derivation.

ORGAN_FREQUENCIES = {
    # Source: Barbara Hero, International Lambdoma Research Institute
    # Base frequency: 256 Hz (scientific C4)
    "Stomach":            110.00,
    "Pancreas":           117.30,
    "Gall Bladder":       164.30,
    "Colon":              176.00,
    "Lungs":              220.00,
    "Intestines":         281.00,
    "Fat Cells":          295.80,
    "Brain":              315.80,
    "Liver":              317.83,
    "Kidneys":            319.88,
    "Blood":              321.90,
    "Muscles":            324.00,
    "Bladder":            352.00,
    "Bone":               418.30,
    "Adrenals & Thyroid": 492.80,
}

ORGAN_FREQUENCIES_META = {
    "researcher": "Barbara Hero",
    "institution": "International Lambdoma Research Institute, Kennebunk, ME",
    "base_frequency": 256,  # Hz, scientific C
    "method": (
        "Sound waves passed through each organ; optimal resonant frequency "
        "calculated using Lambdoma matrix (Pythagorean ratio grid) and "
        "speed of sound. Mathematical derivation, not direct measurement "
        "of organ oscillation."
    ),
    "tradition": "Western alternative / neo-Pythagorean",
    "evidence": "E",
    "note": (
        "Human organs do not have fixed resonant frequencies in the way "
        "that a tuning fork does. Organs are soft tissue with variable "
        "density, blood flow, and boundary conditions. These frequencies "
        "are mathematical constructs, not measured natural resonances."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EXTENDED FREQUENCY COMPENDIUM
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Additional frequencies from various traditions, compiled from the
# healing frequency scale (0-4225 Hz). Cross-references multiple systems.

EXTENDED_FREQUENCIES = {
    # Schumann
    7.83:    {"category": "Schumann",   "claim": "Earth's fundamental EM resonance"},

    # Otto tuning forks
    32:      {"category": "Otto Fork",  "claim": "Immune system stimulation"},
    40:      {"category": "Otto Fork",  "claim": "Nervous system stimulation"},
    64:      {"category": "Brain",      "claim": "Memory and cognition"},

    # Ohm forks
    68.05:   {"category": "Ohm Low",    "claim": "Earth star chakra (Vasundhara); stomach"},
    128:     {"category": "Ohm Mid",    "claim": "Blood circulation, sleep, stress relief; joints; heart chakra"},
    272.2:   {"category": "Ohm High",   "claim": "Soul star chakra (Vyapini)"},
    543.4:   {"category": "Ohm Ultra",  "claim": "Ultra high Ohm"},

    # Angelic frequencies
    111:     {"category": "Angel",      "claim": "Holy frequency; pancreas"},
    222:     {"category": "Angel",      "claim": "Energy balancer"},
    333:     {"category": "Angel",      "claim": "Angelic frequency"},
    444:     {"category": "Angel/Min",  "claim": "Silica mineral frequency"},
    555:     {"category": "Angel",      "claim": "Mental balance improvement"},
    666:     {"category": "Angel",      "claim": "Physical/spiritual balance"},
    777:     {"category": "Angel",      "claim": "Anxiety and nervousness relief"},
    888:     {"category": "Angel",      "claim": "Positive energy and clarity"},
    999:     {"category": "Angel",      "claim": "Higher self frequency"},

    # Tesla frequencies (multiples of 3/6/9)
    324:     {"category": "Tesla",      "claim": "Healing, balance, harmony (3+2+4=9)"},
    639:     {"category": "Tesla/Solf", "claim": "Relationships, communication (6+3+9=18->9)"},
    963:     {"category": "Tesla/Solf", "claim": "Spiritual awakening, divine connection (9+6+3=18->9)"},

    # DNA nucleotide frequencies (unverified)
    528:     {"category": "DNA",        "claim": "Cytosine nucleotide frequency"},
    537.8:   {"category": "DNA",        "claim": "Thymine nucleotide frequency"},
    544.4:   {"category": "DNA",        "claim": "Adenine nucleotide frequency"},
    545.6:   {"category": "DNA",        "claim": "Guanine nucleotide frequency"},

    # Mineral frequencies (from BioAcoustic tradition)
    256:     {"category": "Mineral",    "claim": "Sulphur; cell growth; tissue regeneration"},
    272:     {"category": "Mineral",    "claim": "Selenium and Chlorine"},
    304:     {"category": "Mineral",    "claim": "Potassium"},
    312:     {"category": "Mineral",    "claim": "Platinum"},
    319.88:  {"category": "Mineral",    "claim": "Calcium (also kidneys in Hero system)"},
    336:     {"category": "Mineral",    "claim": "Magnesium"},
    352:     {"category": "Mineral",    "claim": "Sodium, Silver"},
    376:     {"category": "Mineral",    "claim": "Chromium"},
    396:     {"category": "Mineral",    "claim": "Manganese (also Solfeggio UT)"},
    400:     {"category": "Mineral",    "claim": "Iron"},
    418.3:   {"category": "Mineral",    "claim": "Iodine"},
    448:     {"category": "Mineral",    "claim": "Copper"},
    464:     {"category": "Mineral",    "claim": "Phosphorus and Zinc"},

    # High frequencies
    1024:    {"category": "Healing",    "claim": "Energy balancing, immune, pain relief"},
    1152:    {"category": "High",       "claim": "Spiritual enlightenment, transcendence"},
    2172:    {"category": "Angel",      "claim": "Jacob's Ladder"},
    4096:    {"category": "Angel",      "claim": "Pillar of Light"},
    4160:    {"category": "Angel",      "claim": "Stairway to Heaven"},

    # Verdi / alternative concert pitch
    432:     {"category": "Concert",    "claim": "Verdi's A; 'natural' tuning; "
                                                 "C4=256 Hz system; calming properties claimed"},
}

EXTENDED_META = {
    "evidence": "E (vast majority)",
    "note": (
        "Most frequencies in this compendium lack peer-reviewed evidence. "
        "The mineral frequencies, DNA nucleotide frequencies, and angelic "
        "frequencies are practitioner claims without published measurement data. "
        "They are included here for COMPLETENESS of the frequency map, "
        "not as endorsement."
    ),
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# BRAINWAVE BANDS (Reference)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

BRAINWAVE_BANDS = {
    "delta":  {"range_hz": (0.5, 4.0),   "state": "Deep sleep, unconscious repair"},
    "theta":  {"range_hz": (4.0, 8.0),   "state": "Meditation, REM, creativity, shamanic trance"},
    "alpha":  {"range_hz": (8.0, 13.0),  "state": "Relaxed wakefulness, calm focus"},
    "beta":   {"range_hz": (13.0, 30.0), "state": "Active thinking, problem solving"},
    "gamma":  {"range_hz": (30.0, 100.0),"state": "High cognitive processing, binding"},
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CROSS-TRADITION FREQUENCY CONVERGENCE TABLE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# Where multiple traditions converge on the same (or very close)
# frequency, that convergence is notable even if individual claims
# are unverified.

CONVERGENCE_POINTS = [
    {
        "frequency_hz": 136.10,
        "traditions": [
            "Cousto planetary (Earth year)",
            "Indian classical (Sa in some traditions)",
            "Tibetan OM",
        ],
        "note_musical": "C#",
        "significance": "Three independent traditions converge on this frequency",
    },
    {
        "frequency_hz": 7.83,
        "traditions": [
            "Schumann resonance (measured EM phenomenon)",
            "Alpha/Theta brainwave border",
            "Multiple meditation traditions target this state",
        ],
        "significance": "Physically measured + neurologically significant frequency band",
    },
    {
        "frequency_hz_range": (4.0, 4.5),
        "traditions": [
            "Shamanic drumming (cross-cultural: Siberian, Native American, Nordic)",
            "Theta brainwave band",
            "Tibetan singing bowl beat frequencies",
        ],
        "significance": "Repetition rate, not acoustic pitch; cross-cultural convergence",
    },
    {
        "frequency_hz": 256,
        "traditions": [
            "Scientific pitch (C4=256)",
            "Pythagorean tuning fork base",
            "Chinese Gong tone (Earth element)",
            "Lambdoma matrix base (Barbara Hero)",
        ],
        "significance": "Power of 2 (2^8); multiple systems use as reference",
    },
    {
        "frequency_hz": 432,
        "traditions": [
            "Verdi tuning (A=432)",
            "New Age 'natural frequency' claim",
            "Cousto Earth year: A=432.10 Hz when C#=136.10",
        ],
        "significance": "The Cousto connection gives 432 Hz a mathematical (not mystical) basis",
    },
    {
        "frequency_hz": 528,
        "traditions": [
            "Solfeggio (MI)",
            "DNA repair claim (unverified)",
            "Chlorophyll absorption peak claim (partial truth: ~430 nm, not 528 nm)",
        ],
        "significance": "Most commercially promoted healing frequency; least evidence for specific claims",
    },
]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EVIDENCE HIERARCHY SUMMARY
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

EVIDENCE_SUMMARY = {
    "strongest_evidence": [
        "Didgeridoo for sleep apnea (RCT, BMJ 2006, p=0.03)",
        "Shamanic drumming -> theta entrainment (multiple EEG studies)",
        "Schumann resonance exists and overlaps EEG bands (physics + neuroscience)",
        "Indian raga listening -> cortisol reduction, parasympathetic activation",
        "Tibetan bowls produce measurable binaural beats (acoustic physics)",
    ],
    "moderate_evidence": [
        "Chinese FPMT -> measurable effects on stress/insomnia",
        "Sufi dhikr -> alpha wave induction",
        "Pink noise -> slow-wave sleep enhancement (separate from solfeggio)",
        "Binaural beats -> brainwave frequency following",
    ],
    "no_evidence": [
        "528 Hz repairs DNA",
        "Specific Hz values heal specific organs",
        "Solfeggio frequencies derive from Gregorian chant",
        "Mineral frequencies (sulphur=256 Hz, etc.)",
        "Angelic frequencies (111, 222, 333, etc.)",
        "DNA nucleotide frequencies (528, 537.8, 544.4, 545.6 Hz)",
    ],
    "mathematically_valid_but_therapeutically_unproven": [
        "Cousto planetary frequencies (math checks out; healing claims don't)",
        "Lambdoma organ frequencies (mathematical model, not measured resonance)",
        "Pythagorean interval ratios (real math; healing overlay is modern)",
        "Golden ratio / phi frequency relationships",
    ],
}
