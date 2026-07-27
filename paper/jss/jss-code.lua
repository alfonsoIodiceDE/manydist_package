-- Render executable code with the environments supplied by jss.cls.
-- Prompted R source is italicized as CodeInput; printed results remain upright
-- as CodeOutput. Blank source lines are kept without an empty "R>" prompt.

if not quarto.doc.is_format("pdf") then
  return {}
end

local function without_empty_prompts(text)
  local lines = {}

  for line in (text .. "\n"):gmatch("(.-)\n") do
    if line:match("^R>%s*$") then
      table.insert(lines, "")
    else
      table.insert(lines, line)
    end
  end

  return table.concat(lines, "\n")
end

return {
  {
    CodeBlock = function(block)
      local environment = "CodeOutput"

      if block.text:match("^R>") then
        environment = "CodeInput"
      end

      local text = without_empty_prompts(block.text)

      return pandoc.RawBlock(
        "latex",
        "\\begin{" .. environment .. "}\n"
          .. text
          .. "\n\\end{" .. environment .. "}"
      )
    end
  }
}
