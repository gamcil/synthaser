import React, { useMemo } from 'react'

import RuleItem from './Item'


export const RuleList = ({
  domains,
  rules,
  handleAdd,
  handleRemove,
  handleChange,
}) => {

  const domainOptions = useMemo(
    () => domains.map(domain => ({
      label: domain.name,
      value: domain.name
    })), [domains]
  )

  const getItems = () => rules.map(rule => {
    const handleRemoveRule = handleRemove(rule.uuid)
    const handleChangeRule = handleChange(rule.uuid)
    const handleUpdateRule = (label, value) => {
      handleChangeRule({
        target: { name: label, value : value }
      })
    }
    return (
      <RuleItem
        key={rule.uuid}
        name={rule.name}
        domains={rule.domains}
        renames={rule.renames}
        filters={rule.filters}
        evaluator={rule.evaluator}
        allDomains={domains}
        domainOptions={domainOptions}
        handleRemove={handleRemoveRule}
        handleChange={handleChangeRule}
        handleUpdate={handleUpdateRule}
      />
    )
  })

  const items = getItems()

  return (
    <div>
      <div>
        <button type="button" onClick={handleAdd}>Add</button>
      </div>
      <ul>
        {items}
      </ul>
    </div>
  )
}
