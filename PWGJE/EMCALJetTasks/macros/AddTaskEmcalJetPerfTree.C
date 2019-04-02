PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetPerfTree* AddTaskEmcalJetPerfTree(
    const char * suffix = ""
)
{  
  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetPerfTree * task = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetPerfTree::AddTaskEmcalJetPerfTree(suffix);
  return task;
}
