WITH morph_join AS (
    SELECT
    i.individual_code,
    i.individual_sex,
    lh.date as dob,
    p.individual_id,
    p.processing_timestamp,
    p.capture_timestamp,
    m.femoral_abduction_deg,
    m.femoral_adduction_deg,
    m.hip_extension_deg,
    m.hip_external_rotation_deg,
    m.hip_flexion_deg,
    m.hip_internal_rotation_deg
FROM individuals i
    LEFT JOIN individual_life_histories lh
        ON i.individual_id = lh.individual_id
    INNER JOIN processings p
        ON i.individual_id = p.individual_id
    INNER JOIN morphometric_data m
        ON p.processing_id = m.processing_id
WHERE lh.life_history_event = 'birth'
)
SELECT
    individual_code,
    individual_sex,
    EXTRACT(YEAR FROM AGE(processing_timestamp::DATE,dob)) as age,
    EXTRACT(YEAR FROM AGE(processing_timestamp::DATE,dob)) - AVG(EXTRACT(YEAR FROM AGE(processing_timestamp::DATE,dob))) OVER (PARTITION BY individual_code) AS within_age,
    AVG(EXTRACT(YEAR FROM AGE(processing_timestamp::DATE,dob))) OVER (PARTITION BY individual_code) AS between_age,
    femoral_abduction_deg,
    femoral_adduction_deg,
    hip_extension_deg,
    hip_external_rotation_deg,
    hip_flexion_deg,
    hip_internal_rotation_deg,
    dob,
    processing_timestamp,
    capture_timestamp
FROM morph_join
GROUP BY individual_code, individual_sex, femoral_abduction_deg, femoral_adduction_deg, hip_extension_deg, hip_external_rotation_deg, hip_flexion_deg,
         hip_internal_rotation_deg, dob, processing_timestamp, capture_timestamp
ORDER BY individual_code;

