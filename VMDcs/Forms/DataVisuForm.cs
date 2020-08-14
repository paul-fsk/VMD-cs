using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace VMDcs.Forms
{
    public partial class DataVisuForm : Form
    {
        public DataVisuForm()
        {
            InitializeComponent();
        }

        public void setChart(List<Chart> ch)
        {
            this.Controls.Clear();
            int i = 0;
            foreach (var c in ch)
            {
                c.Location = new Point(0, i * 185);
                this.Controls.Add(c);
                i++;
            }

        }
    }
}
