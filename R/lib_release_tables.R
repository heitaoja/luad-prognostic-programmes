# Additional evidence sheets shared by export and verification.
release_evidence_sheets <- function(full) {
  sheets <- c("S16 Driver and smoking" = "03d_driver_smoking_adjustment",
              "S17 Published signatures" = "03e_published_signature_benchmark",
              "S18 IPD and external meta" = "04b_ipd_and_pan_external_synthesis",
              "S19 GEO recalibration" = "05c_recalibration",
              "S20 Signature transport" = "05d_published_signature_transport")
  if (full) sheets <- c(sheets, "S21 CPTAC recalibration" = "05b_rnaseq_recalibration")
  sheets
}
