using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class Phase
{

    public enum Phases : int
    {
        /// <summary>
        /// Subcritical liquid
        /// </summary>
        Liquid = 0,

        /// <summary>
        /// Supercritical (p>pc, T>Tc)
        /// </summary>
        Supercritical,

        /// <summary>
        /// Supercritical gas (p&lt;pc, T&gt;Tc)
        /// </summary>
        SupercriticalGas,

        /// <summary>
        /// Supercritical liquid (p&gt;pc, T&lt;Tc)
        /// </summary>
        SupercriticalLiquid,

        /// <summary>
        /// At critical point
        /// </summary>
        CriticalPoint,

        /// <summary>
        /// Subcritical gas
        /// </summary>
        Gas,

        /// <summary>
        /// Twophase
        /// </summary>
        Twophase,

        /// <summary>
        /// Unknown phase
        /// </summary>
        Unknown,

        NotImposed
    }
}
