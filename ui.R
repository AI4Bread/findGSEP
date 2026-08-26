library(shiny)
library(shinyjs)
library(DT)
library(shinycssloaders)
#library(shinyanimate)
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  #withAnim(),
  tags$head(tags$link(type="text/css", rel="stylesheet", href="header.css"),
            tags$title('findGSEP--Estimating genome size of a species (including polyploids) using k-mer frequencies'),
	     tags$style(HTML(".main-panel {
		background-color: #f5f5f5;
          padding: 10px;
          border: 1px solid #e3e3e3;
          border-radius: 3px;
          box-shadow: 0 0 6px rgba(0, 0, 0, 0.1);
          margin: 30px;
        }"))),
  navbarPage(title="findGSEP Webserver v1.3",id="tabs",position = 'fixed-top', inverse = T,
             tabPanel('Home', icon = icon('home'),
		      fluidRow(
                      column(width=4,
			   img(src = 'findGSEP_logo1205-2.png', height=200,width=300)
			   #imageOutput("UFoldImage"),
			#tags$head(tags$style(
                 	#type="text/css",
                	#"#UFoldImage img {float: center; max-width: 100%; width: 100%; height: auto}"
              		#	))
			),
		      column(width = 8,
			br(),
			br(),
			tags$div(
			tags$h2("Estimating genome size of a species (including polyploids) using k-mer frequencies",style="text-align:center")#,
		        #tags$h2("with Deep Learning",style="text-align:center", class= "grey")
			)),
		      column(width=12,
	    	      tags$div(
	    	        class = "panel panel-primary",
	    	        tags$div(
	    	          class = "panel-heading",
	    	          tags$h4("Usage: ")
	    	        ),
		        tags$style("
 			   .panel-heading {
 			     background-color: #f5f5f5;  /* 设置背景颜色为与well面板相同的颜色 */
 			   }
 			 "),
	    	        tags$div(
	    	          class = "panel-body",
			  tags$p("Upload results from running Jellyfish. Example: ",
              tags$a(href = "javanica.histo", target = "_blank", "javanica_tetraploid.histo"),
              tags$a(href = "O_21mer_new20231012.histo", target = "_blank", "potato_tetraploid.histo"),
		          tags$a(href = "wheat.histo", target = "_blank", "wheat_hexaploid.histo"),
              tags$a(href = "strawberry.histo", target = "_blank", "strawberry_octoploid.histo"),
              tags$a(href = "zhonghuaxun_newgenerate.histo", target = "_blank", "Chinese_sturgeon_octoploid.histo")),

			   tags$h4("Instructions for running Jellyfish:"),
			   tags$ol(
			     tags$li("Download and install jellyfish from: ",
			             tags$a(href = "http://www.genome.umd.edu/jellyfish.html#Release", target = "_blank", "http://www.genome.umd.edu/jellyfish.html#Release")),
			     tags$li("Count k-mers using jellyfish:",
			             tags$p(tags$pre("$ jellyfish count -C -m 21 -t 1 -s 5G *.fastq -o reads.mer")),
			             tags$p("Note you should adjust the memory (-s) and threads (-t) parameter according to your server. This example will use 1 threads and 5GB of RAM. The k-mer length (-m) may need to be scaled if you have low coverage or a high error rate.")),
			     tags$li("Generate the k-mer count histogram file:",
			             tags$p(tags$pre("$ jellyfish histo -h 30000000 -t 10 -o reads.histo reads.mer")),
			             tags$p("Again the thread count (-t) should be scaled according to your server.")),
			     tags$li("Upload reads.histo to findGSEP")
			   ),
			   tags$h4("Instructions for install findGSEP package:"),
			   tags$ol(
                             tags$li("Install devtools package(Under R environment):",
                                     tags$p(tags$pre("$ install.packages(\"devtools\")"))),
                             tags$li("Install directly from github",
                                     tags$p(tags$pre("$ devtools::install_github(\"sperfu/findGSEP\")")))
                           ),
			     tags$li("More instructions can be referred from our ",tags$a(href = "https://github.com/sperfu/findGSEP", target = "_blank", "github repository."))
			   #tags$p("Note: High copy-number DNA such as chloroplasts can confuse the model. Set a max k-mer coverage to avoid this. Default is -1 meaning no filter.")
			 )
			  #markdown::renderMarkdown("instructions.md")
			  #p('2. Enter your output file name (no suffix needed) and click the \'Submit\' button;')
	    	        )
	    	      )
		      ),
		      
                      #br(),
                      #titlePanel('UFold: Fast and Accurate RNA Secondary Structure Prediction with Deep Learning'),
                      #br(),
                      #wellPanel(
                      #  p('This is the k-mer based genome heterozygosity estimation query webserver', class="red"),
                        #p('All licenses in this website are copyrighted by their respective authors. Everything else is released under MIT licenses. '),
                      #  h4('Usage:',class="bold"),
                      #  p('1. Please upload your fasta sequence file or type in fasta sequence in the input box, you can click the \'Show example data\' button below to see an example;'), 
                      #  p('2. Enter your output file name (no suffix needed) and click the \'Submit\' button;'), 
                      #  p('3. Once the job is finished, you can visualize and download your result from the Download panel;'),
                      #  p('* The computation may take some time based on the size of your data. If you have any questions, please feel free to contact us.')
                      #),
                      wellPanel(
                        textAreaInput('your_text', label = "Paste your histo content here:"),
                        selectInput(
                      	'example_sel', 'Select example input', choices = c("javanica_tetraploid","potato_tetraploid","wheat_hexaploid","strawberry_octoploid","Chinese_sturgeon_octoploid"),
                      	 selected = "javanica_tetraploid",
                         selectize = FALSE
                      	),
                        # shinycssloaders::withSpinner(actionButton('example', label = "Use example data")),
                        actionButton('example', label = "Use example data"),	
                        # br(),
                        br(),
                        #textInput('name', label = 'Paste your sequence here:'),
                        fileInput('file_input',label = "Choose file to upload", accept = c("text/csv",".txt",".histo",".hist","text/plain"))
                      ),
                      # wellPanel(
                      #   fluidRow(
                      #     column(numericInput('num', 'decay value', min = 0.001, max=0.1, value = 0.001,width = 150),width = 4),
                      #     column(numericInput('num', 'decay value', min = 0.001, max=0.1, value = 0.001,width = 150),width = 4)
                      #   ),
                      #   radioButtons('extTable', 'Table output format',choices = c('TXT'='txt',"pickle"='pickle'), inline = T)
                      # ),
		      wellPanel(
                        textInput('species_name', label = 'Enter your species name: '),
                        textInput("outfile", label = "Output file name (no extension)", value = "result"),
                        numericInput("n", "Enter ploidy number:", value = 2)
                      ),
                      #passwordInput('passwd', label = 'Password:'),
                      #actionButton('login', label = 'Login'),
                      tags$head(tags$script(src = "message-handler.js")),
                      shinyjs::useShinyjs(),
                      actionButton('submit', label = "Submit", icon = icon('upload')),
                      br(),
                      tableOutput('form_table'),
                      #h5('* MIT License'),
                      br(),
                      br(),
		      #p(id = 'user_num', em("You are current user number: "),style="text-align:center; font-family: times")
		      #p(id = 'user_num', textOutput("numPatientsOutput"),style="text-align:center; font-family: times")
                      # h5('[',a(' Personal Information Collection Statement', href="http://47.104.164.42:3838/sample-apps/ufold/"),'|',
                      #    a(' Github', href="http://47.104.164.42:3838/sample-apps/ufold/"),'|',
                      #    a(' Homepage ', href="http://47.104.164.42:3838/sample-apps/ufold/"),']', class="center"),
             ),

             tabPanel('Download', icon=icon('download'),
             titlePanel('Estimate genome size from sequencing reads using a k-mer-based approach'),
             br(),
             sidebarLayout(
             sidebarPanel(width = 4,
               p('You can download your result from here:'),
               actionButton('draw_plot', label = "Show all data"),
               #downloadButton('ResultDown',label = 'Download Result'),
               #selectInput('seq_name','Choose Seq', c('', 'iris')),
               downloadButton('dVisualize', label = "Download image"),
               #p(id = 'double_click','Please click the Visualize button again if nothing pops up, if you are using Mac OS, you may need to allow pop-up window in the searching tab.',class='red')
	       br(),
	       br(),
               tableOutput('seq_table')
               ),
             
             mainPanel(column(width = 8,offset = 2,
               imageOutput("myImage"),
	       div(
	       style = "margin-top: 50px;",
	       DT::dataTableOutput("table")
		),
               tags$head(tags$style(
                 type="text/css",
                "#myImage img {float: center; width: 100%; height: 100%}",
		"#table {width: 130%;}"
              ))
		)
             )
	     ),
             br(),
             # htmlOutput("inc"),
             # br(),
             #tags$head(tags$style(
             #    type="text/css",
             #   "#panel_varna img {max-width: 100%; width: 100%; max-height: 100%; height: auto}"
             #wellPanel(id='panel_varna',
             #  imageOutput("myImage")
             #),
             #imageOutput("myImage"),
             #wellPanel(id="panel_forna",width = 4,
             #  h4('Acknowledgement', class="bold"),
             #  p('We would like to thank ViennaRNA Web Services Forna for providing the web API to visualize the RNA secondary structure online.'),
             #  h4('Forna citation:',class="bold"),
             #  p('Kerpedjiev P, Hammer S, Hofacker IL (2015). Forna (force-directed RNA): Simple and effective online RNA secondary structure diagrams. Bioinformatics 31(20):3377-9.'),
             #  p('We apologize for any inconvenience caused by the system suspension. If you have any queries, please feel free to contact us.')
             #)
             ),
             #tabPanel('Search', icon=icon('search')),

             tabPanel('Contact', icon = icon('user'),
                      h2('Contact Us'),
                      h3('Contact us if you have any questions about the usage of this website.'),
                      br(),
                      # bachelar
                      h4('Contact information',class='bold'),
                      p('Laiyi Fu & Hequan Sun'),
                      p('Email: fulaiyi@xjtu.edu.cn & hequan.sun@xjtu.edu.cn'),
                      p('Github Page: ', a('Github web', href="https://github.com/sperfu/findGSEP", target="_blank")),
                      #p('To apply: ', a('www.admo.cityu.edu.hk/apply', href='https://www.admo.cityu.edu.hk/apply')),
                      br(),
                      #Msc
                      h4('Cite:',class='bold'),
                      p('Please consider to cite our paper if you use findGSEP in your research.'),
                      p('Paper title:',a('findGSEP: a web application for estimating genome size of polyploid species using k-mer frequencies',href='',target="_blank")),
                      #p('URL: ', a(href='https://www.cityu.edu.hk/sgs/tpg/admission')),
                      p('Corresponding author e-mail: hequan.sun@xjtu.edu.cn'),
                      p('Corresponding author website: ', a('Website', href="https://github.com/HeQSun",target="_blank")),
                      br(),
                      h4('Acknowledgement', class="bold"),
                      p('We would like to thank XiaMeng Wei(Xijing University) for helpful discussion and logo design.'),
                      br(),

                      h4(a('Privacy Policy', href="http://146.56.237.198:3838/findGSEP/"),'|',
                         a('Copyright', href="http://146.56.237.198:3838/findGSP"),'|',
                         a('Disclaimer', href="http://146.56.237.198:3838/findGSEP/")),
                      helpText('The information provided is subject to ongoing review. The University reserves the right to amend the information from time to time.'),
                      helpText('Best Viewed with resolution of 1024 x 768 (or higher) and supports Internet Explorer 8+, Mozilla Firefox 3.5+, Google Chrome and Safari.'),
                      helpText('© 2023 Xi\'an Jiaotong University. All Rights Reserved.')

             )
  )
))
