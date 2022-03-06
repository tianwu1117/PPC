library(shiny)
library(knitr)
library(markdown)
source("RSqrdNorm.R")
source("PowerQuanSs.R")
source("Exp.No.Exp.Var.R")
source("ProbSs.R")
source("samvar.R") #to calculate the phenotypic variance of extreme samples

knit("1.Intro_b.Rmd","1.Intro_b")
#Define UI for application
ui <- navbarPage("Polygenic Power Calculator",
            tabPanel("Introduction", withMathJax(includeMarkdown("1.Intro_b"))
                 ),
      #*******************************************
      navbarMenu("Key GWAS results",
      #*******************************************
            #######################################
            tabPanel("Quantitative phenotype",
            #######################################
            h4("Key GWAS results for quantitative phenotypes"),
                fluidRow(
                  column(4,
                    wellPanel(
                       h4("Key parameters"),
                       numericInput("m.p.quan",  label = "Total number of nearly independent SNP(s) (m)", value = 60000,
                                     min = 10^4, max = 10^7),

                       numericInput("pi0.p.quan",label = "Proportion of null SNPs (\u03C0\U2080)", value = 0.99,
                                     min = 0.9, max = 0.9999),

                       numericInput("h2.p.quan", label = "SNP heritability (h\U00B2)",value = 0.4,
                                     min = 0.01, max = 0.9),

                       numericInput("n.p.quan",  label = "Training sample size (n)",value = 10^5,
                                     min = 10^4,max = 10^8),

                       numericInput("alpha.p.quan",label = "Significance level (\u03B1)",value = 5*10^-8,
                                     min = 10^-8,max = 0.5),

                       numericInput("p.given.p.quan",label = "Probability of detecting at least X significant SNP(s)  ", value = 0.9,
                                     min=0.01,  max = 1),

                       numericInput("x.p.quan","The number of  X significant SNP(s)",
                                     value = 10, min = 1),
                       # p(strong("Please finalize your parameters before plotting because it willtake a while:")),
                       # actionButton("QQplot.quan", HTML("Click this button and Submit <br/>
                       #                                   to plot the expected QQ plot "), style = 'width: 240px' ),
                       #
                       # actionButton("TDRplot.quan",HTML("Click this button and Submit <br/>
                       #               to plot the epected local TDR plot"), style = 'width: 240px'),
                       submitButton("Submit")
                    )
                  ),
                  column(8,
                        h3(p("Results:")),
                        htmlOutput("WarningMsg"),
                        htmlOutput("Exp.No.Sig"),    #expected no. of independent significant SNP(s)
                        htmlOutput("SE.no.sig.SNP"), #standard error of the no. of independent significant SNP(s)
                        htmlOutput("exp.var.exp"),      #expected variance explained by the independent significant SNP(s)
                        htmlOutput("SE.var.exp"),       #standard error of the variance explained by significant SNP(s)
                        tags$br(),
                        htmlOutput("Exp.No.True.Sig"),   #expected no. of independent significant SNP(s)
                        htmlOutput("SE.no.true.sig.SNP"),#standard error of the expected no. of TRUE independent SNP(s)
                        htmlOutput("Exp.true.var.exp"),     #expected variance explained by the TRUE independent significant SNP(s)
                        htmlOutput("SE.true.var.exp"),      #standard error of the variance explained by the TRUE significant SNPs
                        tags$br(),
                        htmlOutput("Power.Quan.1.ss"),     #sample size
                        htmlOutput("Power.Quan.2.prob"),   #probability
                        #plotOutput("Power.Quan.5.QQ"),     #expected QQ plot
                        #plotOutput("Power.Quan.5.TDR")     #expected local TDR plot
                        tags$style(type="text/css",
                                   ".shiny-output-error { visibility: hidden; }",
                                   ".shiny-output-error:before { visibility: hidden; }"
                        )
                  )
                )
            ),
            #######################################
            tabPanel("Extreme sample selection",
            #######################################
            h4("Key GWAS results under extreme sample selection design"),
              fluidRow(
                column(3,
                       "Please first calculate the phenotypic variance of extreme sample",
                       wellPanel(
                         numericInput("TL.p", "Lower threshold for extreme sample selection",
                                      value = -1.96),
                         numericInput("TU.p", "Upper threshold for extreme sample selection",
                                      value =  1.96),
                         numericInput("ratio.p.extr","Proportion of samples below the lower threshold",
                                      value = 0.5),
                         submitButton("Submit")
                       ),
                         htmlOutput("Power.extr.sampvar") #sample variance of the extreme samples
                      ),
                column(4,
                 wellPanel(
                   numericInput("m.p.extr", label = "Total number of nearly independent SNP(s) (m)", value = 60000,
                                 min = 10^4,max = 10^7),

                   numericInput("pi0.p.extr", label = "Proportion of null SNPs (\u03C0\U2080)", value = 0.99,
                                 min = 0.9, max = 0.9999),

                   numericInput("h2.p.extr", label = "SNP heritability (h\U00B2)",value = 0.4,
                                 min = 0.01,max = 0.9),

                   numericInput("n.p.extr", label = "Training sample size (n)", value = 10000,
                                min = 1000, max = 10^7),

                   numericInput("alpha.p.extr",label = "Significance level (\u03B1)", value = 5*10^-8,
                                 min = 10^-10,max = 0.5),

                   numericInput("varY.p.extr","Phenotypic variance of extreme sample",
                                value = 5.58, min = 1),

                   numericInput("p.given.p.extr",label = "Probability of detecting at least X significant SNP(s)  ", value = 0.9,
                                 min=0.01,  max = 1),

                   numericInput("x.p.extr","The number of  X significant SNP(s)", value = 10,
                                 min = 1),
                   # p(strong("Please finalize your parameters before plotting because it will take a while:")),
                   # actionButton("QQplot.extr",HTML("Click this button and Submit <br/>
                   #                to plot the expected QQ plot "), style = 'width: 240px'),
                   # actionButton("TDRplot.extr",HTML("Click this button and Submit <br/>
                   #                   to plot the epected local TDR plot"), style = 'width: 240px'),
                   submitButton("Submit")
                         )
                    ),
                column(5,
                  h3(p("Results:")),
                  htmlOutput("Extr.Exp.No.Sig"),    #expected no. of independent significant SNP(s)
                  htmlOutput("Extr.SE.no.sig.SNP"), #standard error of the no. of independent significant SNP(s)
                  htmlOutput("Extr.exp.var.exp"),      #expected variance explained by the independent significant SNP(s)
                  htmlOutput("Extr.SE.var.exp"),       #standard error of the variance explained by significant SNP(s)
                  tags$br(),
                  htmlOutput("Extr.Exp.No.True.Sig"),   #expected no. of independent significant SNP(s)
                  htmlOutput("Extr.SE.no.true.sig.SNP"),#standard error of the expected no. of TRUE independent SNP(s)
                  htmlOutput("Extr.Exp.true.var.exp"),     #expected variance explained by the TRUE independent significant SNP(s)
                  htmlOutput("Extr.SE.true.var.exp"),      #standard error of the variance explained by the TRUE significant SNPs
                  tags$br(),
                  htmlOutput("Power.extr.1.ss"),     #sample size
                  htmlOutput("Power.extr.2.prob")   #probability
                  # plotOutput("Power.extr.5.QQ"),
               # plotOutput("Power.extr.5.TDR")
                       )
              )
            ), #tabPanel("Extreme sample selection") done

            #######################################
            tabPanel("Binary phenotype",
            #######################################
            h4("Key GWAS results for case-control study"),
                  sidebarPanel(
                    numericInput("m.p.bi",label = "Total number of nearly independent SNP(s) (m)", value = 60000,
                                 min = 10^4,max = 10^7),

                    numericInput("pi0.p.bi","Proportion of null SNPs (\u03C0\U2080)",
                                 value = 0.99),

                    numericInput("h2.p.bi",label = "SNP heritability (h\U00B2)", value = 0.99,
                                 min = 0.9, max = 0.9999),

                    numericInput("n.p.bi", label = "Training sample size", value = 10000,
                                 min = 1000, max = 10^7),

                    numericInput("alpha.p.bi",label = "Significance level (\u03B1)",value = 5*10^-8,
                                 min = 10^-10,max = 0.5),

                    numericInput("K.p.bi","Disease prevalence in population (K)", 0.1,
                                 min = 0.001, max = 0.5),

                    numericInput("w.p.bi","Proportion of cases in training sample (w)", value = 0.5,
                                 min = 0.01,  max = 100),

                    numericInput("p.given.p.bi","Probability of detecting at least X significant SNP(s)", value = 0.9,
                                 min=0.01, max = 1),

                    numericInput("x.p.bi","The number of  X significant SNP(s)", value = 10,
                                  min = 1),

                    # p(strong("Please finalize your parameters before plotting because it will take a while:")),
                    # actionButton("QQplot.bi",HTML("Click this button and Submit <br/>
                    #                                      to plot the expected QQ plot "), style = 'width: 240px'),
                    # br(),
                    # actionButton("TDRplot.bi",HTML("Click this button and Submit <br/>
                    #                  to plot the epected local TDR plot"), style = 'width: 240px'),
                    # br(),
                    submitButton("Submit")
                                      ),
                      mainPanel(
                        h3(p("Results:")),
                        htmlOutput("Bi.Exp.No.Sig"),    #expected no. of independent significant SNP(s)
                        htmlOutput("Bi.SE.no.sig.SNP"), #standard error of the no. of independent significant SNP(s)
                        htmlOutput("Bi.exp.var.exp"),      #expected variance explained by the independent significant SNP(s)
                        htmlOutput("Bi.SE.var.exp"),       #standard error of the variance explained by significant SNP(s)
                        tags$br(),
                        htmlOutput("Bi.Exp.No.True.Sig"),   #expected no. of independent significant SNP(s)
                        htmlOutput("Bi.SE.no.true.sig.SNP"),#standard error of the expected no. of TRUE independent SNP(s)
                        htmlOutput("Bi.Exp.true.var.exp"),     #expected variance explained by the TRUE independent significant SNP(s)
                        htmlOutput("Bi.SE.true.var.exp"),      #standard error of the variance explained by the TRUE significant SNPs
                        tags$br(),
                        htmlOutput("Power.bi.1.ss"),     #sample size
                        htmlOutput("Power.bi.2.prob")   #probability
                              )
                           )   #tabPanel("Binary phenotype") done
                  ), #navbarMenu("Key GWAS results") done

      ############################################
      #*******************************************
      navbarMenu("Prediction accuracy",
      #*******************************************
            #######################################
            tabPanel("Quantitative phenotype",
            #######################################
            h4("Prediction accuracy of quantitative phenotype"),
            #Sidebar layout with a input and output definitions
                  sidebarLayout(
                # Inputs: select values of parameters
                     sidebarPanel(
                          numericInput("m","Total number of nearly independent SNP(s) (m)", value = 60000,
                                        min = 10^4, max = 10^7),

                          numericInput("pi0","Proportion of null SNPs (\u03C0\U2080)", value = 0.99,
                                      min = 0.8, max = 0.999),

                          numericInput("h2","SNP heritability (h\U00B2)",value = 0.4,
                                      min = 0.1, max = 0.99),

                          numericInput("n","Training sample size (n)", value = 10^5,
                                      min = 1000, max = 10^6),
                          submitButton("Submit")
                     ),
                # Output: Show R^2 of OLSE and posterior expectation
                     mainPanel(
                          htmlOutput("DescpRand"),
                          tags$br(),
                          tableOutput("RSqrdTNoparaRandSample") #R2 calculated by method(s) without tuning parameter method for random population sample
                              )
                                                    )
                                ),  #tabPanel("Quantitative phenotype") done
            #######################################
            tabPanel("Extreme sample selection",
            #######################################
            h4("Prediction accuracy of extreme sample selection"),
                  sidebarLayout(
                     sidebarPanel(
                          numericInput("m.extr",label = "Total number of nearly independent SNP(s) (m)", value = 60000,
                                        min = 10^4, max = 10^7),

                          numericInput("pi0.extr","Proportion of null SNPs (\u03C0\U2080)", value = 0.99,
                                      min = 0.8, max = 0.999),

                          numericInput("h2.extr","SNP heritability (h\U00B2)", value = 0.4,
                                      min = 0.1, max = 0.99),

                          numericInput("n.extr","Training sample size (n)", value = 10^5,
                                      min = 1000, max = 10^6),

                          numericInput("varY.extr","Phenotypic variance of extreme sample",
                                               value = 1.5, min = 1, max = 10),

                          submitButton("Submit")
                                ),
                     mainPanel(
                        htmlOutput("samvarmsg"),
                        tags$br(),
                        tableOutput(outputId ="RSqrdTNoparaExtreSample"))
                               )
                     ), #tabPanel("Extreme sample selection") done
            #######################################
            tabPanel("Binary phenotype",
            #######################################
            h4("Prediction accuracy of binary phenotype"),
                  sidebarLayout(
                     sidebarPanel(
                          numericInput("m.bi",label = "Total number of nearly independent SNP(s) (m)", value = 60000,
                                        min = 10^4, max = 10^7),

                          numericInput("pi0.bi","Proportion of null SNPs (\u03C0\U2080):", value = 0.99,
                                      min = 0.8, max = 0.999),

                          numericInput("h2.bi","SNP heritability (h\U00B2)", value = 0.4,
                                      min = 0.1, max = 0.99),

                          numericInput("n.bi","Training sample size (n)", value = 10^5,
                                      min = 1000, max = 10^8),

                          numericInput("K.bi","Disease prevalence in population (K)", 0.1,
                                       min = 0.001, max = 0.5),

                          numericInput("w.bi","Proportion of cases in training sample (w)", value = 0.5,
                                       min = 0.01,  max = 100),

                         submitButton("Submit")
                                     ),
                     mainPanel(
                          htmlOutput("DescpBin"),
                          tags$br(),
                          tableOutput(outputId ="RSqrdTNoparaLinear")
                              )
                                  )
                        )  #tabPanel("Binary phenotype") done
                      ),#navbarMenu("Prediction accuracy") done

            #######################################
            tabPanel("Reference",
                     p("Should you have any questions of this power calculator, please refer to the manuscript"),
                     p("Polygenic Power Calculator: Statistical Power and Polygenic Prediction Accuracy of Genome-wide Association Studies of Complex Traits."),
                     tags$br(),
                     p("Contact person: Talia Wu, u3006094 at connect.hku.hk")
                   )
            ) #navbarPage done
            #######################################

server <- function(input, output){

  output$DescpRand <- renderText({
   "Prediction accuracy of polygenic score is measured by R\U00B2, which is the squared correlation between the true and estimated polygenic scores."
  })
##############################################################
#1. Results of Prediction accuracy of PGS
##############################################################
#1.1 Quantitative phenotype

MetdPool1 <- c("OLSE","Posterior mean")

 output$RSqrdTNoparaRandSample <- renderTable({
    RSqrd.pool <- vector(length = 2, mode = "numeric")
    RSqrd.pool <- RSqrdNorm(input$h2, input$m,input$n,input$pi0)
    df <- data.frame(MetdPool1,RSqrd.pool)
    names(df) <- c("Methods of constructing PGS", "R\U00B2")
    df
    }, digits = 2, width = "400px")


#1.2 Extreme case selction
 #message about sample variance calculation for extreme sample design
 output$samvarmsg <- renderText({
   "If you have not calculated the variance of your extreme sample, please go to the 'Key GWAS results'-
   'extreme sample selection' page to calculate it first."
 })

 output$RSqrdTNoparaExtreSample <- renderTable({

     RSqrd.pool.extr <- RSqrdNorm(input$h2.extr, input$m.extr,input$n.extr*(input$varY.extr)^2,input$pi0.extr)
     df <- data.frame(MetdPool1,RSqrd.pool.extr)
     names(df) <- c("Method", "R\U00B2")
     df}, digits = 6, width = "255px")

 #1.3 Binary phenotype
 output$DescpBin <- renderText({
    "The R\U00B2 is the squared correlation coefficient between the estimated and true liability of disease."
 })

 output$RSqrdTNoparaLinear <- renderTable({
     fac <- (input$w.bi*(1-input$w.bi)*dnorm(qnorm(input$K.bi,lower.tail = F))^2)/(input$K.bi*(1-input$K.bi))^2
     RSqrd.pool.bi <- RSqrdNorm(input$h2.bi, input$m.bi,input$n.bi*fac,input$pi0.bi)
     df <- data.frame(MetdPool1,RSqrd.pool.bi)
     names(df) <- c("Method", "R\U00B2")
     df}, digits = 6,  width = "255px")

##############################################################
#2. Key GWAS results
##############################################################

    ##########################################################
    #check the validity of input
    ##########################################################
    ###m###
    observeEvent(c(input$m.p.quan,input$m.p.extr,input$m.p.bi,input$m,input$m.extr,input$m.bi), {
        if(input$m.p.quan >= 10^8 | input$m.p.quan <= 0 | !is.numeric(input$m.p.quan)|
           input$m.p.extr >= 10^8 | input$m.p.extr <= 0 | !is.numeric(input$m.p.extr)|
           input$m.p.bi >= 10^8 | input$m.p.bi <= 0 | !is.numeric(input$m.p.bi)      |
           input$m >= 10^8 | input$m <= 0 | !is.numeric(input$m)                     |
           input$m.extr >= 10^8 | input$m.extr <= 0 | !is.numeric(input$m.extr)      |
           input$m.bi >= 10^8 | input$m.bi <= 0 | !is.numeric(input$m.bi))
           {
        showModal(modalDialog(
          title = "Warning: Invalid input",
          "m must be a number between 0 and 1e8."))
           }
        })
    ###h2###
    observeEvent(c(input$h2.p.quan,input$h2.p.extr,input$h2.p.bi,input$h2,input$h2.extr,input$h2.bi), {
    if(input$h2.p.quan >= 1 | input$h2.p.quan <= 0 | !is.numeric(input$h2.p.quan)|
       input$h2.p.extr >= 1 | input$h2.p.extr <= 0 | !is.numeric(input$h2.p.extr)|
       input$h2.p.bi >= 1 | input$h2.p.bi <= 0 | !is.numeric(input$h2.p.bi)|
       input$h2 >= 1 | input$h2 <= 0 | !is.numeric(input$h2)|
       input$h2.extr >= 1 | input$h2.extr <= 0 | !is.numeric(input$h2.extr)|
       input$h2.bi >= 1 | input$h2.bi <= 0 | !is.numeric(input$h2.bi))
    {
       showModal(modalDialog(
          title = "Warning: Invalid input",
          "h\U00B2 must be a number between 0 and 1."))
    }
 })
    ###pi0###
    observeEvent(c(input$pi0.p.quan,input$pi0.p.extr,input$pi0.p.bi,input$pi0,input$pi0.extr,input$pi0.bi), {
    if(input$pi0.p.quan >= 10^8 | input$pi0.p.quan <= 0 | !is.numeric(input$pi0.p.quan)|
       input$pi0.p.extr >= 10^8 | input$pi0.p.extr <= 0 | !is.numeric(input$pi0.p.extr)|
       input$pi0.p.bi >= 10^8 | input$pi0.p.bi <= 0 | !is.numeric(input$pi0.p.bi)|
       input$pi0 >= 10^8 | input$pi0 <= 0 | !is.numeric(input$pi0)|
       input$pi0.extr >= 10^8 | input$pi0.extr <= 0 | !is.numeric(input$pi0.extr)|
       input$pi0.bi >= 10^8 | input$pi0.bi <= 0 | !is.numeric(input$pi0.bi))
    {
       showModal(modalDialog(
          title = "Warning: Invalid input",
          "\u03C0\U2080 must be a number between 0 and 1."))
    }
 })
    ###n###
    observeEvent(c(input$n.p.quan,input$n.p.extr,input$n.p.bi,input$n,input$n.extr,input$n.bi), {
       if(input$n.p.quan >= 10^8 | input$n.p.quan <= 0 | !is.numeric(input$n.p.quan)|
          input$n.p.extr >= 10^8 | input$n.p.extr <= 0 | !is.numeric(input$n.p.extr)|
          input$n.p.bi >= 10^8 | input$n.p.bi <= 0 | !is.numeric(input$n.p.bi)|
          input$n >= 10^8 | input$n <= 0 | !is.numeric(input$n)|
          input$n.extr >= 10^8 | input$n.extr <= 0 | !is.numeric(input$n.extr)|
          input$n.bi >= 10^8 | input$n.bi <= 0 | !is.numeric(input$n.bi)
          )
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "n must be a number between 0 and 1e8."))
       }
    })
    ###alpha###
    observeEvent(c(input$alpha.p.quan,input$alpha.p.extr,input$alpha.p.bi), {
       if(input$alpha.p.quan >= 10^8 | input$alpha.p.quan <= 0 | !is.numeric(input$alpha.p.quan)|
          input$alpha.p.extr >= 10^8 | input$alpha.p.extr <= 0 | !is.numeric(input$alpha.p.extr)|
          input$alpha.p.bi >= 10^8 | input$alpha.p.bi <= 0 | !is.numeric(input$alpha.p.bi))
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "\u03B1 must be a number between 0 and 1."))
       }
    })
    ###given probability###
    observeEvent(c(input$p.given.p.quan,input$p.given.p.extr,input$p.given.p.bi), {
       if(input$p.given.p.quan >= 1 | input$p.given.p.quan <= 0 | !is.numeric(input$p.given.p.quan)|
          input$p.given.p.extr >= 1 | input$p.given.p.extr <= 0 | !is.numeric(input$p.given.p.extr)|
          input$p.given.p.bi >= 1 | input$p.given.p.bi <= 0 | !is.numeric(input$p.given.p.bi))
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "The probability must be a number between 0 and 1."))
       }
     })
    ###Number of SNPs###
    observeEvent(c(input$x.p.quan,input$x.p.extr,input$x.p.bi), {
       if(input$x.p.quan <= 0 | !is.numeric(input$x.p.quan) | (round(input$x.p.quan) != input$x.p.quan)|
          input$x.p.extr <= 0 | !is.numeric(input$x.p.extr) | (round(input$x.p.extr) != input$x.p.extr)|
          input$x.p.bi <= 0 | !is.numeric(input$x.p.bi) | (round(input$x.p.bi) != input$x.p.bi))
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "The number of significant SNPs must be an integer between 0 and the total number of non-null SNPs."))
       }
    })
    ###Phenotypic variance of extreme sample###
    observeEvent(c(input$varY.extr, input$varY.p.extr),{
       if(input$varY.extr <1   | !is.numeric(input$varY.extr) |
          input$varY.p.extr <1 | !is.numeric(input$varY.p.extr))
          {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "The phenotypic variace under the extreme sample design should be greater than 1."))
       }
    })
    ###Disease prevalence###
    observeEvent(c(input$K.bi, input$K.p.bi),{
       if(input$K.bi <= 0 |input$K.bi >= 1  |!is.numeric(input$K.bi)|
          input$K.p.bi <= 0 | input$K.p.bi >= 1 |!is.numeric(input$K.p.bi))
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "Diease prevalence must be a number between 0 and 1."))
       }
    })
    ###Proportion of cases###
    observeEvent(c(input$w.bi, input$w.p.bi),{
       if(input$w.bi <= 0 |input$w.bi >= 1  |!is.numeric(input$w.bi)|
          input$w.p.bi <= 0 | input$w.p.bi >= 1 |!is.numeric(input$w.p.bi))
       {
          showModal(modalDialog(
             title = "Warning: Invalid input",
             "The proportion of cases must be a number between 0 and 1."))
       }
    })
    ###Ratio of extreme small samples
    observeEvent(input$ratio.p.extr,{
      if(input$ratio.p.extr <=0 | input$ratio.p.extr>=1 | !is.numeric(input$ratio.p.extr)){
        showModal(modalDialog(
          title = "Warning: Invalid input",
          "The proportion must be a number between 0 and 1."))
      }
    })
    ###Bound of extreme samples
    observeEvent(c(input$TL.p,input$TU.p),{
      if(input$TL.p >= input$TU.p |!is.numeric(input$TL.p) |!is.numeric(input$TU.p)){
        showModal(modalDialog(
          title = "Warning: Invalid input",
          "The upper bound of extremely small samples must be a number less than the lower bound of extremely large samples."))
      }
    })


 output$Exp.No.Sig <- renderText({
   paste0("The expected number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[1],"</b></font>", ".")
 })

 output$SE.no.sig.SNP <- renderText({
   paste0("The standard error of the number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
         Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[2],"</b></font>", ".")
 })

 output$exp.var.exp <- renderText({
   paste0("The expected variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[3],"</b></font>",".")
 })

 output$SE.var.exp <- renderText({
   paste0("The standard error of the variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[4],"</b></font>",".")

 })

 output$Exp.No.True.Sig <- renderText({
   paste0("The expected number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[5],"</b></font>",".")
 })

 output$SE.no.true.sig.SNP <- renderText({
   paste0("The standard error of the number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[6],"</b></font>",".")
 })

 output$Exp.true.var.exp <- renderText({
   paste0("The expected variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[7],"</b></font>",".")
 })

 output$SE.true.var.exp <- renderText({
   paste0("The standard error of the variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan,input$alpha.p.quan)[8],"</b></font>",".")
 })

 output$Power.Quan.1.ss <- renderText({

    paste0("Given ", round(input$p.given.p.quan * 100,2),"% probability to detect at least " , round(input$x.p.quan) ,
           " truly associated SNP(s), the desirable sample size is " ,"<font color=\"#0000FF\"><b>",
           formatC(PowerQuanSs(input$p.given.p.quan, input$x.p.quan, input$alpha.p.quan, input$h2.p.quan, input$m.p.quan,input$pi0.p.quan),
                   format = "f", big.mark = ",", digits = 0),"</b></font>",".")
 })

 output$Power.Quan.2.prob <- renderText({
    paste0("Given sample size is ", formatC(input$n.p.quan, format = "f", big.mark = ",", digits = 0),
           " , the probability that you could detect at least ", round(input$x.p.quan) ," truly associated SNP(s) is " ,"<font color=\"#0000FF\"><b>",
           ProbExpNo(input$alpha.p.quan, input$x.p.quan, input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan),"</b></font>",".")
 })
    # }
#})



 # QQ.quan <- eventReactive(input$QQplot.quan, {
 #   withProgress(message = 'Plotting the expected QQ plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.quan/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan)$p.val
 # })

 # TDR.quan <- eventReactive(input$TDRplot.quan, {
 #   withProgress(message = 'Plotting the expected TDR plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.quan/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.quan, input$m.p.quan, input$n.p.quan, input$pi0.p.quan)$tdr
 # })
 #
 # output$Power.Quan.5.QQ <- renderPlot({
 #
 #   # par(mfrow = c(1,2))
 #   expected <- (1:input$m.p.quan)/(1+input$m.p.quan)
 #   plot(-log10(expected),-log10(sort(QQ.quan())),
 #         xlab = "-log10(expected p-value under NULL)",
 #         ylab = "-log10(expected observed p-value under model)",
 #         main = "Expected QQ plot")
 #   abline(0,1)
 # })
 #
 # output$Power.Quan.5.TDR <- renderPlot({
 #   expected <- (1:input$m.p.quan)/(1+input$m.p.quan)
 #   plot(sort(-log10(expected)),TDR.quan(),
 #        xlab = "-log10(expected p-value under NULL)",
 #        ylab = "local TDR",
 #        main = "local TDR vs expected p-values under NULL")
 # })

 output$Power.extr.sampvar <- renderText({
   paste0("The variance of your sample is ", "<font color=\"#0000FF\"><b>",
   round(samvar(input$TL.p,input$TU.p,input$ratio.p.extr),2),"</b></font>",". You could proceed to key GWAS results calculation by putting
          this number in the 'Phenotypic variance of extreme sample' box and inputting other paremeters.")
 })

 output$Extr.Exp.No.Sig <- renderText({
   paste0("Given sample size is ", formatC(input$n.p.extr, format="f", big.mark=",", digits=0) ,
          " , the expected number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[1],"</b></font>",".")
 })

 output$Extr.SE.no.sig.SNP <- renderText({
   paste0("The standard error of the number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[2],"</b></font>",".")
 })

 output$Extr.exp.var.exp <- renderText({
   paste0("The expected variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[3],"</b></font>",".")
 })

 output$Extr.SE.var.exp <- renderText({
   paste0("The standard error of the variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[4],"</b></font>",".")
 })

 output$Extr.Exp.No.True.Sig <- renderText({
   paste0("The expected number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[5],"</b></font>",".")
 })

 output$Extr.SE.no.true.sig.SNP <- renderText({
   paste0("The standard error of the number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
         Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[6],"</b></font>",".")
 })

 output$Extr.Exp.true.var.exp <- renderText({
   paste0("The expected variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
         Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[7],"</b></font>",".")
 })

 output$Extr.SE.true.var.exp <- renderText({
   paste0("The standard error of the variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr,input$alpha.p.extr)[8],
          "</b></font>",".")
 })

 output$Power.extr.1.ss <- renderText({
    paste0("Given ", round(input$p.given.p.extr * 100,2),"% probability to detect at least ", round(input$x.p.extr) , " truly associated SNP(s), the desirable sample size is " ,
           "<font color=\"#0000FF\"><b>",
           formatC(PowerQuanSs(input$p.given.p.extr, input$x.p.extr, input$alpha.p.extr, input$h2.p.extr, input$m.p.extr,input$pi0.p.extr)/(input$varY.p.extr)^2,
                   format="f", big.mark=",", digits=0),"</b></font>",".")
 })

 output$Power.extr.2.prob <- renderText({
    paste0("Given sample size is ", formatC(input$n.p.extr, format="f", big.mark=",", digits=0),
           " , the probability that you could detect at least ", input$x.p.extr ," truly associated SNP(s) is " ,"<font color=\"#0000FF\"><b>",
           ProbExpNo(input$alpha.p.extr, input$x.p.extr, input$h2.p.extr, input$m.p.extr, input$n.p.extr*(input$varY.p.extr)^2, input$pi0.p.extr),"</b></font>",".")
 })

 # QQ.extr <- eventReactive(input$QQplot.extr, {
 #   withProgress(message = 'Plotting the expected QQ plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.extr/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.extr, input$m.p.extr, input$n.p.extr, input$pi0.p.extr)$p.val
 # })
 #
 # TDR.extr <- eventReactive(input$TDRplot.extr, {
 #   withProgress(message = 'Plotting the expected TDR plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.extr/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.extr, input$m.p.extr, input$n.p.extr, input$pi0.p.extr)$tdr
 # })

 # output$Power.extr.5.QQ <- renderPlot({
 #
 #   # par(mfrow = c(1,2))
 #   expected <- (1:input$m.p.extr)/(1+input$m.p.extr)
 #   plot( -log10(expected),-log10(sort(QQ.extr())),
 #         xlab = "-log10(expected p-value under NULL)",
 #         ylab = "-log10(expected observed p-value under model)",
 #         main = "Expected QQ plot")
 #   abline(0,1)
 # })
 #
 # output$Power.extr.5.TDR <- renderPlot({
 #   expected <- (1:input$m.p.extr)/(1+input$m.p.extr)
 #   plot(sort(-log10(expected)),TDR.extr(),
 #        xlab = "-log10(expected p-value under NULL)",
 #        ylab = "local TDR",
 #        main = "local TDR vs expected p-values under NULL")
 # })

 output$Bi.Exp.No.Sig <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("Given sample size is ", formatC(input$n.p.bi, format="f", big.mark=",", digits=0),
          " and the proportion of cases is ", input$n.w.bi ," , the expected number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[1],"</b></font>",".")
 })

 output$Bi.SE.no.sig.SNP <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The standard error of the number of significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
         Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[2],"</b></font>",".")
 })

 output$Bi.exp.var.exp <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The expected variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[3],"</b></font>",".")
 })

 output$Bi.SE.var.exp <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The standard error of the variance explained by the significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[4],"</b></font>",".")
 })

 output$Bi.Exp.No.True.Sig <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The expected number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          formatC(Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[5],
                  format="f", big.mark=",", digits=0),"</b></font>",
          ".")
 })

 output$Bi.SE.no.true.sig.SNP <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The standard error of the number of TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          formatC(Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[6],
                  format="f", big.mark=",", digits=0),"</b></font>",
          ".")
 })

 output$Bi.Exp.true.var.exp <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The expected variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[7],"</b></font>",".")
 })

 output$Bi.SE.true.var.exp <- renderText({
   fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
   paste0("The standard error of the variance explained by TRUE significant SNP(s) is " ,"<font color=\"#0000FF\"><b>",
          Exp.No.Exp.Var(input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi,input$alpha.p.bi)[8],"</b></font>",".")
 })

 output$Power.bi.1.ss <- renderText({
    fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
    #n real * fac = effective sample size n*
    paste0("Given ", round(input$p.given.p.bi * 100,2),"% probability to detect at least ",round(input$x.p.bi)," truly associated SNP(s), the desirable sample size is " ,
           "<font color=\"#0000FF\"><b>",
           formatC(PowerQuanSs(input$p.given.p.bi, input$x.p.bi, input$alpha.p.bi, input$h2.p.bi, input$m.p.bi,input$pi0.p.bi)/fac,2 ,
                   format="f", big.mark=",", digits=0),"</b></font>",".")
 })

 output$Power.bi.2.prob <- renderText({
    fac <- (dnorm(qnorm(input$K.p.bi, lower.tail = F))^2*input$w.p.bi*(1-input$w.p.bi))/(input$K.p.bi*(1-input$K.p.bi))^2
    paste0("Given sample size is ", formatC(input$n.p.bi, format="f", big.mark=",", digits=0) ,
    " , the probability that you could detect at least ",input$x.p.bi , " truly associated SNP(s) is " ,"<font color=\"#0000FF\"><b>",
           ProbExpNo(input$alpha.p.bi, input$x.p.bi, input$h2.p.bi, input$m.p.bi, input$n.p.bi*fac, input$pi0.p.bi),"</b></font>",".")
 })
 # QQ.bi <- eventReactive(input$QQplot.bi, {
 #   withProgress(message = 'Plotting the expected QQ plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.bi/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.bi, input$m.p.bi, input$n.p.bi, input$pi0.p.bi)$p.val
 # })
 #
 # TDR.bi <- eventReactive(input$TDRplot.bi, {
 #   withProgress(message = 'Plotting the expected TDR plot',
 #                detail = 'This may take a while...', value = 0, {
 #                  for (i in 1:10) {
 #                    incProgress(1/10)
 #                    Sys.sleep(input$m.p.bi/50000)
 #                  }
 #                }
 #   )
 #   PvalTDR(input$h2.p.bi, input$m.p.bi, input$n.p.bi, input$pi0.p.bi)$tdr
 # })

 # output$Power.bi.5.QQ <- renderPlot({
 #
 #   # par(mfrow = c(1,2))
 #   expected <- (1:input$m.p.bi)/(1+input$m.p.bi)
 #   plot( -log10(expected),-log10(sort(QQ.bi())),
 #         xlab = "-log10(expected p-value under NULL)",
 #         ylab = "-log10(expected observed p-value under model)",
 #         main = "Expected QQ plot")
 #   abline(0,1)
 # })
 #
 # output$Power.bi.5.TDR <- renderPlot({
 #   expected <- (1:input$m.p.bi)/(1+input$m.p.bi)
 #   plot(sort(-log10(expected)),TDR.bi(),
 #        xlab = "-log10(expected p-value under NULL)",
 #        ylab = "local TDR",
 #        main = "local TDR vs expected p-values under NULL")
 # })

}

shinyApp(ui = ui, server = server)
