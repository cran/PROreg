\name{EDpro}
\alias{EDpro}
\alias{EDpro}
\title{
Eating Disorders patient-reported outcome data.
}
\description{
A prospective study of patients diagnosed with an eating disorder who were followed up for two years in the Eating Disorders Outpatient Clinic of the Psychiatric Service at Galdakao Hospital in Bizkaia, Basque Country.
}

\usage{
data(EDpro)}


\details{
All consecutive patients seen at the Eating Disorders Outpatient Clinic of the Psychiatric Service at Galdakao Hospital between October 1996 and October 1997 were eligible for the study. Patients were included if they were ambulatory; had received a diagnosis of anorexia nervosa or bulimia nervosa based on criteria in the Diagnostic and Statistical Manual of Mental Disorders, 4th edition (DSM-IV); and were between the ages of 14 and 65 years.

Throughout the 2-year study period, each of the 193 patient received a psychopharmacologic and psychotherapeutic treatment program consisting of cognitive-behavioural treatment; nutritional orientation and counselling; psychoeducation; motivationaltherapy; social skills training; and therapy to modify distorted perception of body image. These interventions were adjusted to each patient's needs by a multidisciplinary team.

At most, 3 measurements were obtained from each of the 193 patients: at baseline, at one year, and at two years from the start of the study. Clinical measurements and Helath-Related Quality of Life (HQRoL) data were recorded at each visit using the Short Form-36 (SF-36) and the Hospital Anxiety and Depression (HAD) Scale helath questionnaires.

More information about the study can be found in Padierna et al. (2002) and Arostegui et al. (2007).}

\format{
The \code{EDpro} data frame has 525 rows and 18 columns. The variables included in the data frame are the following:

-id: Patient identifier.

-visit: Visit number.

-bmi: Body Mass Index.

age: Age.

-duration: Duration of the illness in years.

-diagnosis: Diagnosis of the eating disorder (1: Anorexia Nervosa, 2: Bulimic Anorexia, 3: Bulimia Nervosa/Binge Eating).

-severity: Severity of the disease, using the Eating Attitudes Test (1: Mild, 2: Moderate, 3: Severe).

-anxiety: Level of anxiety using the Hospital Anxiety and Depression Scale (0: No, 1: Yes).

-depression: Level of depression using the Hospital Anxiety and Depression Scale (0: No, 1: Yes).

-pf: Physical Function dimension measured by the SF-36 questionnaire.

-rp: Role Physical dimension measured by the SF-36 questionnaire.

-bp: Body Pain dimension measured by the SF-36 questionnaire.

-gh: General Health dimension measured by the SF-36 questionnaire.

-vt: Vitality dimension measured by the SF-36 questionnaire.

-sf: Social Functioning dimension measured by the SF-36 questionnaire.

-re: Role Emotional dimension measured by the SF-36 questionnaire.

-mh: Mental Health dimension measured by the SF-36 questionnaire.

-time: Years elapsed since the beginning of the study.

}

\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2007): Analysis of short-form-36 (SF-36): The
beta-binomial distribution approach, \emph{Statistics in Medicine}, \bold{26}, 1318-1342.

Padierna A., Quintana J.M., Arostegui I., Gonzalez N. & Horcajo M.J.(2002): Changes in health related quality of life among patients treated for eating disorders, \emph{Quality of Life Research}, \bold{11}, 545-552.
}
