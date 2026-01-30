tags$head(
  tags$style(HTML("
    /* New Shiny progress modal */

    .shiny-notification {
      position: fixed !important;
      top: 50% !important;
      left: 50% !important;
      transform: translate(-50%, -50%) !important;
      width: 600px;
      max-width: 90%;
      text-align: center;
      padding: 30px;
      font-size: 18px;
    }

    .shiny-notification-content {
      text-align: center;
    }

    .progress {
      height: 30px;
      margin-top: 20px;
    }

    .progress-bar {
      font-size: 16px;
      line-height: 30px;
    }
  "))
)
