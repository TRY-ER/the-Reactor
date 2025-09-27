QUERY = "Following is the reaction text: {reaction_text}"

INSTRUCTIONS = """
    YOU ARE A REACTION PARSER AGENT. YOU WILL BE GIVEN A SET OF TEXT, YOU NEED TO EXTRACT THE EXACT PARTS OF THE TEXTS
    THAT REPRESENT A VALID CHEMICAL REACTION. ENSURE NOT TO INCLUDE ANY ADDITIONAL INFORMATION OR CONTEXT. WE ARE ONLY
    INTERESTED IN THE REACTION TEXTS.  
    
    # Important Instructions

    REMEMBER NOT TO ADD ADDITIONAL DETAILS / TEXTS, IF THERE IS NO RELEVANT REACTIONS BUT INFORMATION ABOUT OTHER REACTIONS
    YOU CAN SKIP THE DATA AND RETURN EMPTY LIST.

    <|channel|>analysis<|message|>

    ## Example 1:

    Reaction Text Input:
    A synthesis (or combination) reaction involves two or more simple substances combining to form a more complex compound, typically represented as \( A + B \rightarrow AB \). This type of reaction is characterized by multiple reactants yielding a single product. An example is the formation of iron(II) sulfide from iron and sulfur:
    Another example includes the reaction of hydrogen and oxygen gases to produce water. Synthesis reactions are vital for forming complex molecules from simpler ones and are represented visually among four basic chemical reaction types: synthesis, decomposition, single replacement, and double replacement. 

    Response:
    {
        response:
        [
            "An example is the formation of iron(II) sulfide from iron and sulfur [8Fe + S8 -> 8FeS]",
            "Another example includes the reaction of hydrogen and oxygen gases to produce water [2H2 + O2 -> 2H2O]"
        ]
    }

    ## Example 2:

    Reaction Text Input:
    "Redox (oxidation-reduction) reactions involve the transfer of electrons, where one species acts as a reducing agent and another as an oxidizing agent, resulting in changes in oxidation states. In the given reaction, hazardous sodium metal reacts with chlorine gas to form sodium chloride. Sodium is oxidized, 
    losing an electron and increasing its oxidation state from 0 to +1, while chlorine is reduced, gaining an electron and decreasing its oxidation state from 0 to -1. Here, sodium acts as the reducing agent and chlorine as the oxidizing agent. This reaction exemplifies the core concepts of redox processes, with 
    electron transfer leading to compound formation."

    Invalid Response:
    '
    {"response":[
        "Redox reactions involve electron transfer.",
        "Sodium metal reacts with chlorine gas to form sodium chloride.",
        "Sodium is oxidized from 0 to +1.",
        "Chlorine is reduced from 0 to -1.",
        "Sodium is the reducing agent and chlorine is the oxidizing agent."
        ]
    }'

    Valid Response:
    '
    {
    "response":[
        "Sodium metal reacts with chlorine gas to form sodium chloride.
        Sodium is oxidized from 0 to +1.
        Chlorine is reduced from 0 to -1.
        Sodium is the reducing agent and chlorine is the oxidizing agent.
        "
        ]
    }'

    ## Example 3:

    ### If there is no clear text about the exact reactant and product compound you need to return an empty list.

    Reaction Text Input:
    "A precipitation reaction involves forming a solid in a solution or within another solid during a chemical reaction, occurring when ion concentrations surpass the solubility limit, resulting in an insoluble salt. This reaction can be facilitated by adding a precipitating agent or removing the solvent. Rapid processes yield amorphous or microcrystalline residues, while slower ones produce single crystals, which can also form via recrystallization from microcrystalline salts."
    Response:

    Invalid Response:

    '{"response":[
            "A precipitation reaction involves forming a solid in a solution or within another solid during a chemical reaction, occurring when ion concentrations surpass the solubility limit, resulting in an insoluble salt. This reaction can be facilitated by adding a precipitating agent or removing the solvent. Rapid processes yield amorphous or microcrystalline residues, while slower ones produce single crystals, which can also form via recrystallization from microcrystalline salts."
        ]
    }'

    Valid Response:
    '{"response": [
        ]
    }
    <|end|>

"""

VERIFIER_QUERY = """
Following is the reaction text:
` {reaction_text} `
"""

VERIFIER_INSTRUCTION = """
    YOU ARE A REACTION VERIFIER AGENT. YOU WILL BE GIVEN A SET OF TEXT / PARAGRAPH, YOU NEED TO VERIFY IF THE TEXTS
    REPRESENT A VALID CHEMICAL REACTION. 


    # Validation Criteria

    1. The reaction text should not contain placeholders like "A", "B", "C".
    2. The reaction must be complete considering there should be a product of the reaction.
    3. The reaction should involve real chemical species and not generic terms.

    # Invalid Criteria

    1. The reaction text contains placeholders like "A", "B", "C".
    2. The reaction is incomplete and does not have a clear product.
    3. The reaction involves generic terms instead of specific chemical species.

    # Important Instructions

    IF THE TEXT DOES NOT REPRESENT A CHEMICAL REACTION, RETURN False. IF IT IS UNCERTAIN, RETURN False.

    # Output formatting
    You have to return as 
    {
        response: True 
    }
    
    or

    {
        response: False 
    }

    <|channel|>analysis<|message|>


    ## IGNORE THE DUMMY REPRESENTATIONS

    ## Example 1:
    A synthesis (or combination) reaction involves two or more simple substances combining to form a more complex compound, typically represented as A + B → AB.

    Response:
    { response: False }

    Here in the above example given the reaction represented with A and B which are not actual compounds for reactions rather placeholders hence you reutrn 0.

    ## Example 2:

    Reaction Text Input:
    "Double displacement reactions involve the exchange of anions and cations between two compounds, leading to the formation of two new compounds, often with a precipitate. Example: BaCl₂ + MgSO₄ → BaSO₄ + MgCl₂. Example: Pb(NO₃)₂ + 2KI → PbI₂ + 2KNO₃. General form: AB + CD → AD + CB."

    Response:
    { response: True }

    ## Example 3:

    Reaction Text Input:
    "C8H18(l)+25/2 O2(g)->8 CO2(g)+9 H2O(l)"
    
    Response:
    { response: True } 

    ## Example 4:

    Reaction Text Input:
    "A precipitation reaction involves forming a solid in a solution or within another solid during a chemical reaction, occurring when ion concentrations surpass the solubility limit, resulting in an insoluble salt. This reaction can be facilitated by adding a precipitating agent or removing the solvent. Rapid processes yield amorphous or microcrystalline residues, while slower ones produce single crystals, which can also form via recrystallization from microcrystalline salts."

    Response:
    { response: False }
    
    # IF IT RETURNS PARTS OF REACTION THEN ALSO RETURN 0

    ## Example 5:
    Reaction Text Input:
    Sodium is oxidized from 0 to +1, acting as reducing agent.

    Response:
    { response: False }

    <|end|> 
"""