function ShowProgress(app, title, perc)
    app.StatusPanel.Title = strcat('Status: ', title);
    numBars = uint8(perc/10);

    switch(numBars)
        case 0
            app.perc0.Visible = true;
        case 1
            app.perc1.Visible = true;
        case 2
            app.perc2.Visible = true;
        case 3
            app.perc3.Visible = true;
        case 4
            app.perc4.Visible = true;
        case 5
            app.perc5.Visible = true;
        case 6
            app.perc6.Visible = true;
        case 7
            app.perc7.Visible = true;
        case 8
            app.perc8.Visible = true;
        case 9
            app.perc9.Visible = true;
        case 10
            app.perc0.Visible = false;
            app.perc1.Visible = false;
            app.perc2.Visible = false;
            app.perc3.Visible = false;
            app.perc4.Visible = false;
            app.perc5.Visible = false;
            app.perc6.Visible = false;
            app.perc7.Visible = false;
            app.perc8.Visible = false;
            app.perc9.Visible = false;
            app.StatusPanel.Title = 'Status: ';
    end
    drawnow
end