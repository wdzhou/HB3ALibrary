using System;
using System.Collections.Generic;
using System.Text;

namespace DAS.UBMatrix
{
    /// <summary>
    /// Unit Cell storage and operation
    /// </summary>
    public class UnitCell
    {
        #region private properties
        double[] axiles;
        double[] angles;
        double volume;
        /// <summary>
        /// 1 for direct/realspace cell; -1 for reciprocal cell
        /// </summary>
        int celltype;
        #endregion

        #region Constructors
        /// <summary>
        /// construct a unit cell
        /// </summary>
        public UnitCell(){
            axiles = new double[3];
            angles = new double[3];
            volume = -1;
            celltype = 0;
            return;
        }
        #endregion

        #region public properties
        /// <summary>
        /// Type of cell, real space or reciprocal space
        /// </summary>
        public int Type
        {
            get { return celltype; }
            set { celltype = value; }
        }
        /// <summary>
        /// a
        /// </summary>
        public double a{
            get {return axiles[0];}
            set {axiles[0] = value;}
        }
        /// <summary>
        /// b
        /// </summary>
        public double b{
            get {return axiles[1];}
            set {axiles[1] = value;}
        }
        /// <summary>
        /// c
        /// </summary>
        public double c{
            get {return axiles[2];}
            set {axiles[2] = value;}
        }
        /// <summary>
        /// alpha
        /// </summary>
        public double alpha{
            get {return angles[0];}
            set {angles[0] = value;}
        }
        /// <summary>
        /// beta
        /// </summary>
        public double beta{
            get {return angles[1];}
            set {angles[1] = value;}
        }
        /// <summary>
        /// gamma
        /// </summary>
        public double gamma{
            get {return angles[2];}
            set {angles[2] = value;}
        }
        /// <summary>
        /// vector as [a, b, c]
        /// </summary>
        public double[] Axiles
        {
            get { return axiles; }
        }
        /// <summary>
        /// vector as [alpha, beta, gamma]
        /// </summary>
        public double[] Angles
        {
            get { return angles; }
        }
        /// <summary>
        /// unit cell volume
        /// </summary>
        public double Volume
        {
            set { volume = value; }
            get { return volume; }
        }
            
        #endregion

        #region Public methods

        /// <summary>
        /// Calculate the volume of the unit cell
        /// (usually called after all axiles and angles are set)
        /// </summary>
        public void CalcualteVolume(out string logmessage)
        {
            // 1. Calculate sin(angles) and cos(angles)
            double[] sangles = new double[3];
            for (int i = 0; i < 3; i++)      
                sangles[i] = DasMath.SinD(angles[i]);            

            double[] cangles = new double[3];
            for (int i = 0; i < 3; i++)
                cangles[i] = Math.Sqrt(1.0 - sangles[i] * sangles[i]);

            // 2. Calculate volume
            // volume=directcell->a*directcell->b*directcell->c
            //        *sqrt(1.0-calpha*calpha-cbeta*cbeta-cgamma*cgamma+2.0*calpha*cbeta*cgamma);
            double tempv1 = 0;
            for (int i = 0; i < 3; i++)
            {
                tempv1 += cangles[i] * cangles[i];
            }
            double tempv2 = 1;
            for (int i = 0; i < 3; i++)
            {
                tempv2 *= cangles[i];
            }
            volume = a * b * c * Math.Sqrt(1.0 - tempv1 + 2.0 * tempv2);

            // Console.WriteLine("Volume = {0} In CalVolume()", volume);

            logmessage = String.Format("Calculate Unit Cell Volume: ");
            for (int i = 0; i < 3; i++)
                logmessage += String.Format("{0} ", axiles[i]);
            for (int i = 0; i < 3; i++)
                logmessage += String.Format("{0} ", angles[i]);
            logmessage += String.Format("\nVolume = {0}", volume);
           
           
            return;
        }

        /// <summary>
        /// from the directcell parameters this calculates the recprical cell and the reciprical volume.
        /// </summary>
        /// <param name="rcell"></param>
        public void CalculateReciprocalUnitCell(UnitCell rcell)
        {
            // 1. Check whether volume has been calculated
            string logmsg; 
            if (volume <= 0)
                CalcualteVolume(out logmsg);

            double[] sangles = new double[3];
            for (int i = 0; i < 3; i++)
                sangles[i] = DasMath.SinD(angles[i]);

            double[] cangles = new double[3];
            for (int i = 0; i < 3; i++)
                cangles[i] = DasMath.CosD(angles[i]);       
            // cangles[i] = Math.Sqrt(1.0 - sangles[i] * sangles[i]);

            // 3. Construct the reciprocal unitcell 
            // Console.WriteLine("volume = {0} In CalculateReciprocal...()", volume);
            rcell.a = (b * c * sangles[0]) / volume;
            rcell.b = (a * c * sangles[1]) / volume;
            rcell.c = (a * b * sangles[2]) / volume;
            rcell.alpha = DasMath.AcosD((cangles[1] * cangles[2] - cangles[0]) / (sangles[1] * sangles[2]));
            rcell.beta  = DasMath.AcosD((cangles[2] * cangles[0] - cangles[1]) / (sangles[0] * sangles[2]));
            rcell.gamma = DasMath.AcosD((cangles[0] * cangles[1] - cangles[2]) / (sangles[1] * sangles[0]));

            if (volume != 0.0)
                rcell.volume = 1.0 / volume;
            else
                rcell.volume = 1.0E20;                

            return;
        }

        /// <summary>
        /// Calculate direct unit cell from G matrix or
        /// Calculate reciprocal unit cell from G^-1 matrix
        /// </summary>
        /// <param name="mG">3x3 matrix as G natrux</param>
        public void CalculateFromG(Matrix mG)
        {
            // 1. a, b, c
            for (int i = 0; i < 3; i++)
                axiles[i] = Math.Sqrt(mG.M[i][i]);

            // 2. alpha, beta, gamma
            alpha = DasMath.AcosD(mG.M[1][2] / (b * c));
            beta  = DasMath.AcosD(mG.M[0][2] / (a * c));
            gamma = DasMath.AcosD(mG.M[0][1] / (a * b));

            // 3. Volume
            string logmsg;
            CalcualteVolume(out logmsg);

            return;
        }

        /// <summary>
        /// Calculate metric matrix (G) from unit cell
        /// </summary>
        /// <source>SUBROUTINE METR(CELL,GI,VOL2)</source>
        /// <param name="mG">[out] Matrix G</param>
        public void CalculateG(Matrix mG)
        {
            // 1. contruct matrix G from input cell
            for (int i = 0; i < 3; i++)
                mG.M[i][i] = axiles[i] * axiles[i];
#if false
            mG.M[0][0] = cell.a * cell.a;
            mG.M[1][1] = cell.b * cell.b;
            mG.M[2][2] = cell.c * cell.c;
#endif
            mG.M[0][1] = a * b * DasMath.CosD(gamma);
            mG.M[1][0] = mG.M[0][1];
            mG.M[0][2] = a * c * DasMath.CosD(beta);
            mG.M[2][0] = mG.M[0][2];
            mG.M[1][2] = b * c * DasMath.CosD(alpha);
            mG.M[2][1] = mG.M[1][2];

            return;
        }

        /// <summary>
        /// Copy the current lattice parameters and type to the target unit cell
        /// </summary>
        /// <param name="targetcell"></param>
        public void CopyTo(UnitCell targetcell)
        {
            for (int i = 0; i < 3; i++)
            {
                targetcell.axiles[i] = axiles[i];
                targetcell.angles[i] = angles[i];
            }
            targetcell.celltype = celltype;
            targetcell.volume = volume;
            return;
        }


        /// <summary>
        /// Print the information of this unit cell
        /// </summary>
        public void Print()
        {
            string ostr = ToString();
            Console.WriteLine(ostr);

            return;
        }

        /// <summary>
        /// Customerized format output
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            string logmsg;
            string ostr = "";
            ostr += String.Format("a, b, c = [{0}, {1}, {2}]  ",
                a.ToString("0.00000"), b.ToString("0.00000"), c.ToString("0.00000"));
            ostr += String.Format("alpha, beta, gamma = [{0}, {1}, {2}]  ",
                alpha.ToString("0.00000"), beta.ToString("0.00000"), gamma.ToString("0.00000"));
            if (volume <= 0)
                CalcualteVolume(out logmsg);
            ostr += String.Format("Volume = {0}", volume);

            return ostr;
        }

        #endregion

        #region static public method

        /// <summary>
        /// SUBROUTINE METR(CELL,GI,VOL2)
        /// Obtain reciprocal matrical matrix and square of volume
        /// </summary>
        /// <param name="cell"></param>
        /// <param name="mGI"></param>
        /// <param name="vol2"></param>
        static public void GenMetricMatrix(UnitCell cell, Matrix mGI, out double vol2)
        {
            Matrix mG = new Matrix(3);

            // 1. contruct matrix G from input cell
            mG.M[0][0] = cell.a * cell.a;
            mG.M[1][1] = cell.b * cell.b;
            mG.M[2][2] = cell.c * cell.c;
            mG.M[0][1] = cell.a * cell.b * DasMath.CosD(cell.gamma);
            mG.M[1][0] = mG.M[0][1];
            mG.M[0][2] = cell.a * cell.c * DasMath.CosD(cell.beta);
            mG.M[2][0] = mG.M[0][2];
            mG.M[1][2] = cell.b * cell.c * DasMath.CosD(cell.alpha);
            mG.M[2][1] = mG.M[1][2];

            // 2. Obtain reciprical metrical matrix and the volume squared.
            int errcode = mG.InvertMatrix(mGI, out vol2);

            return;
        }

        #endregion

    }
}